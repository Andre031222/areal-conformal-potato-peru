import math
import pathlib
import re
import urllib.request

import rasterio
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

RAIZ = pathlib.Path(__file__).resolve().parents[1] / "data"
OESTE, SUR, ESTE, NORTE = -81.5, -18.6, -68.5, 0.5


SUELO = RAIZ / "soilgrids"
SUELO.mkdir(exist_ok=True)
for variable in ["phh2o", "soc", "clay", "sand"]:
    for capa in ["0-5cm", "5-15cm", "15-30cm"]:
        destino = SUELO / f"{variable}_{capa}_mean_peru.tif"
        if destino.exists():
            continue
        vrt = f"/vsicurl/https://files.isric.org/soilgrids/latest/data/{variable}/{variable}_{capa}_mean.vrt"
        with rasterio.Env(GDAL_HTTP_MAX_RETRY=4, GDAL_HTTP_RETRY_DELAY=5):
            with rasterio.open(vrt) as fuente:
                caja = transform_bounds("EPSG:4326", fuente.crs, OESTE, SUR, ESTE, NORTE, densify_pts=50)
                ventana = from_bounds(*caja, transform=fuente.transform).round_offsets().round_lengths()
                datos = fuente.read(1, window=ventana)
                perfil = fuente.profile | {
                    "driver": "GTiff", "height": datos.shape[0], "width": datos.shape[1],
                    "transform": fuente.window_transform(ventana), "compress": "deflate", "tiled": True,
                }
                with rasterio.open(destino, "w", **perfil) as salida:
                    salida.write(datos, 1)
        print("suelo", destino.name, datos.shape)


DEM = RAIZ / "dem_glo90"
DEM.mkdir(exist_ok=True)
URL_DEM = ("https://copernicus-dem-90m.s3.amazonaws.com/"
           "Copernicus_DSM_COG_30_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
           "Copernicus_DSM_COG_30_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif")
for lat in range(math.floor(SUR), math.ceil(NORTE)):
    for lon in range(math.floor(OESTE), math.ceil(ESTE)):
        ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
        url = URL_DEM.format(ns=ns, lat=abs(lat), ew=ew, lon=abs(lon))
        destino = DEM / url.rsplit("/", 1)[1]
        if destino.exists():
            continue
        try:
            urllib.request.urlretrieve(url, destino)
        except Exception:
            pass
print("dem", len(list(DEM.glob("*.tif"))), "teselas")


AGRO = RAIZ / "superficie_agricola"
AGRO.mkdir(exist_ok=True)
destino = AGRO / "superficie_agricola_2024.zip"
if not destino.exists():
    identificador = "11XHc0X8fQb1u-XcnIPqPQmDPTp9nPUyU"
    abridor = urllib.request.build_opener(urllib.request.HTTPCookieProcessor())
    url = f"https://drive.usercontent.google.com/download?id={identificador}&export=download&confirm=t"
    with abridor.open(url) as respuesta, open(destino, "wb") as salida:
        salida.write(respuesta.read())
print("superficie agrícola", destino.stat().st_size, "bytes")
