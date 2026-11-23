import pathlib
import time
from concurrent.futures import ThreadPoolExecutor

import geopandas as gpd
import numpy as np
import pandas as pd
import planetary_computer
import pystac_client
import rasterio
from rasterio.features import rasterize
from rasterio.warp import Resampling, reproject

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
SALIDA = DATOS / "modis"
PARTES = SALIDA / "partes"
CACHE = SALIDA / "rejillas"
for carpeta in (SALIDA, PARTES, CACHE):
    carpeta.mkdir(parents=True, exist_ok=True)

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLECCION = "modis-13Q1-061"
INICIO, FIN = "2000-09-01", "2025-09-30"

distritos = gpd.read_file(DATOS / "limites" / "DISTRITO.gpkg")[["ubigeo", "geometry"]]
distritos["id"] = np.arange(1, len(distritos) + 1, dtype=np.int32)
UBIGEO = np.concatenate([[""], distritos["ubigeo"].astype(str).str.zfill(6).to_numpy()])

MASCARA = DATOS / "superficie_agricola" / "mascara_agricola_2024.tif"


def catalogo():
    return pystac_client.Client.open(STAC, modifier=planetary_computer.sign_inplace)


def rejilla(item):
    tesela = f"h{item.properties['modis:horizontal-tile']:02d}v{item.properties['modis:vertical-tile']:02d}"
    ruta = CACHE / f"{tesela}.npz"
    if ruta.exists():
        guardado = np.load(ruta)
        return tesela, guardado["ids"], guardado["frac"]
    with rasterio.open(item.assets["250m_16_days_NDVI"].href, overview_level=0) as r:
        forma, transformada, crs = r.shape, r.transform, r.crs
    zona = distritos.to_crs(crs)
    ids = rasterize(zip(zona.geometry, zona["id"]), out_shape=forma, transform=transformada,
                    fill=0, dtype="int32")

    frac = np.zeros(forma, dtype="float32")
    with rasterio.open(MASCARA) as fuente:
        reproject(fuente.read(1).astype("float32"), frac, src_transform=fuente.transform,
                  src_crs=fuente.crs, dst_transform=transformada, dst_crs=crs,
                  resampling=Resampling.average, src_nodata=None, dst_nodata=0)
    np.savez_compressed(ruta, ids=ids, frac=frac)
    return tesela, ids, frac


def procesar(item, ids, frac):
    destino = PARTES / f"{item.id}.parquet"
    if destino.exists():
        return 0
    for intento in range(5):
        try:
            ruta_ndvi = planetary_computer.sign(item.assets["250m_16_days_NDVI"].href.split("?")[0])
            ruta_calidad = planetary_computer.sign(item.assets["250m_16_days_pixel_reliability"].href.split("?")[0])
            with rasterio.open(ruta_ndvi, overview_level=0) as r:
                ndvi = r.read(1)
            with rasterio.open(ruta_calidad, overview_level=0) as r:
                calidad = r.read(1)
            break
        except Exception as error:
            print("reintento", item.id, error, flush=True)
            time.sleep(15 * (intento + 1))
    else:
        print("FALLÓ", item.id, flush=True)
        return 0
    valido = (ids > 0) & (frac > 0.01) & (ndvi >= -2000) & (ndvi <= 10000) & (calidad >= 0) & (calidad <= 1)
    clave = ids[valido]
    peso = frac[valido]
    suma = np.bincount(clave, weights=peso * ndvi[valido] * 1e-4, minlength=len(UBIGEO))
    total = np.bincount(clave, weights=peso, minlength=len(UBIGEO))
    cuenta = np.bincount(clave, minlength=len(UBIGEO))
    hay = np.nonzero(total)[0]
    pd.DataFrame({
        "ubigeo": UBIGEO[hay],
        "fecha": pd.to_datetime(item.id.split(".")[1][1:], format="%Y%j"),
        "tesela": item.id.split(".")[2],
        "suma": suma[hay], "peso": total[hay], "n": cuenta[hay],
    }).to_parquet(destino, index=False)
    return len(hay)


def main():
    cat = catalogo()

    muestra = list(cat.search(collections=[COLECCION], bbox=list(distritos.total_bounds),
                              datetime="2024-01-01/2024-01-16").items())
    teselas = {}
    for item in muestra:
        tesela, ids, frac = rejilla(item)
        if (ids > 0).any() and (frac > 0).any():
            teselas[tesela] = (ids, frac)
    print("teselas con tierra agrícola del Perú:", sorted(teselas), flush=True)

    for tesela, (ids, frac) in sorted(teselas.items()):
        h, v = int(tesela[1:3]), int(tesela[4:6])
        items = list(cat.search(collections=[COLECCION], datetime=f"{INICIO}/{FIN}",
                                query={"modis:horizontal-tile": {"eq": h},
                                       "modis:vertical-tile": {"eq": v}}).items())
        unicos = {}
        for it in sorted(items, key=lambda x: x.id):
            unicos[it.id.split(".")[1]] = it
        items = list(unicos.values())
        with ThreadPoolExecutor(max_workers=8) as grupo:
            hechos = list(grupo.map(lambda it: procesar(it, ids, frac), items))
        print(tesela, len(items), "compuestos;", sum(1 for x in hechos if x), "nuevos", flush=True)

    elegidos = {}
    for parte in sorted(PARTES.glob("*.parquet")):
        clave = tuple(parte.name.split(".")[1:3])
        elegidos[clave] = parte
    partes = pd.concat(pd.read_parquet(p) for p in elegidos.values())
    junto = partes.groupby(["ubigeo", "fecha"], as_index=False)[["suma", "peso", "n"]].sum()
    junto["ndvi"] = junto["suma"] / junto["peso"]
    junto = junto[(junto["peso"] >= 0.5) & junto["ndvi"].between(-0.2, 1.0)]
    junto[["ubigeo", "fecha", "ndvi", "peso", "n"]].to_parquet(SALIDA / "ndvi_16dias_distrito.parquet", index=False)
    print("filas", len(junto), "distritos", junto["ubigeo"].nunique(),
          "fechas", junto["fecha"].min().date(), junto["fecha"].max().date(), flush=True)


if __name__ == "__main__":
    main()
