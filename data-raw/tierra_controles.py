"""Suelo, relieve, accesibilidad y riego por distrito 2012, sobre los píxeles agrícolas.

- SoilGrids 2.0 (250 m, Homolosine): pH, carbono orgánico, arcilla y arena, promedio 0–30 cm
  ponderado por el espesor de cada capa (5, 10 y 15 cm). Unidades: pH, g/kg y %.
- Copernicus GLO-90: altitud (media, p10, p90) y pendiente en grados, calculada en cada tesela.
- Malaria Atlas Project 2015: minutos a la ciudad más cercana (mediana del distrito).
- IV CENAGRO 2012: fracción de la superficie sembrada bajo riego.

Salida: datos/tierra_controles_distrito.parquet.
"""
import pathlib
import re

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
CAPAS = {"0-5cm": 5, "5-15cm": 10, "15-30cm": 15}
SUELO = {"phh2o": ("ph", 0.1), "soc": ("carbono_org_gkg", 0.1), "clay": ("arcilla_pct", 0.1), "sand": ("arena_pct", 0.1)}

puntos = pd.read_parquet(DATOS / "puntos_agricolas.parquet")
lon, lat = puntos["lon"].to_numpy(), puntos["lat"].to_numpy()
salida = pd.DataFrame({"ubigeo": puntos["ubigeo"]})


def muestrear(ruta, x, y):
    with rasterio.open(ruta) as r:
        banda = r.read(1)
        f, c = rasterio.transform.rowcol(r.transform, x, y)
        f, c = np.asarray(f), np.asarray(c)
        dentro = (f >= 0) & (f < r.height) & (c >= 0) & (c < r.width)
        valor = np.full(len(x), np.nan)
        valor[dentro] = banda[f[dentro], c[dentro]]
        if r.nodata is not None:
            valor[valor == r.nodata] = np.nan
    return valor


with rasterio.open(DATOS / "soilgrids" / "phh2o_0-5cm_mean_peru.tif") as r:
    homolosine = Transformer.from_crs("EPSG:4326", r.crs, always_xy=True)
hx, hy = homolosine.transform(lon, lat)
for prop, (nombre, escala) in SUELO.items():
    suma = np.zeros(len(lon))
    peso = np.zeros(len(lon))
    for capa, espesor in CAPAS.items():
        v = muestrear(DATOS / "soilgrids" / f"{prop}_{capa}_mean_peru.tif", hx, hy) * escala
        ok = np.isfinite(v)
        suma[ok] += v[ok] * espesor
        peso[ok] += espesor
    salida[nombre] = np.where(peso > 0, suma / np.maximum(peso, 1), np.nan)
    print(nombre, "faltantes:", int(np.isnan(salida[nombre]).sum()))

altitud = np.full(len(lon), np.nan)
pendiente = np.full(len(lon), np.nan)
for ruta in sorted((DATOS / "dem_glo90").glob("*.tif")):
    ns, la, ew, lo = re.search(r"_([NS])(\d+)_00_([EW])(\d+)_00_", ruta.name).groups()
    sur = int(la) * (1 if ns == "N" else -1)
    oeste = int(lo) * (1 if ew == "E" else -1)
    dentro = np.nonzero((lat >= sur) & (lat < sur + 1) & (lon >= oeste) & (lon < oeste + 1))[0]
    if not len(dentro):
        continue
    with rasterio.open(ruta) as r:
        z = r.read(1).astype("float64")
        dy = abs(r.transform.e) * 110_574
        lat_f = r.transform.f + (np.arange(r.height) + 0.5) * r.transform.e
        dx = (r.transform.a * 111_320 * np.cos(np.radians(lat_f)))[:, None]
        gy, gx = np.gradient(z)
        pend = np.degrees(np.arctan(np.hypot(gx / dx, gy / dy)))
        f, c = rasterio.transform.rowcol(r.transform, lon[dentro], lat[dentro])
        f = np.clip(np.asarray(f), 0, r.height - 1)
        c = np.clip(np.asarray(c), 0, r.width - 1)
        altitud[dentro] = z[f, c]
        pendiente[dentro] = pend[f, c]
salida["altitud_m"] = altitud
salida["pendiente_grados"] = pendiente
print("altitud faltante:", int(np.isnan(altitud).sum()))

salida["viaje_min"] = muestrear(DATOS / "accesibilidad" / "tiempo_viaje_ciudades_2015_peru.tif", lon, lat)
salida.loc[salida["viaje_min"] < 0, "viaje_min"] = np.nan
salida["area_ha"] = (0.00125 * 111_320) ** 2 * np.cos(np.radians(lat)) / 10_000
salida["agricola"] = puntos["base"].eq("agricola").to_numpy()

g = salida.groupby("ubigeo")
tabla = g[[v[0] for v in SUELO.values()] + ["pendiente_grados"]].mean()
tabla["altitud_media_m"] = g["altitud_m"].mean()
tabla["altitud_p10_m"] = g["altitud_m"].quantile(0.10)
tabla["altitud_p90_m"] = g["altitud_m"].quantile(0.90)
tabla["viaje_min"] = g["viaje_min"].median()
tabla["ha_mascara_agricola"] = salida[salida["agricola"]].groupby("ubigeo")["area_ha"].sum()
tabla["ha_mascara_agricola"] = tabla["ha_mascara_agricola"].fillna(0.0)
tabla["base_agregacion"] = np.where(tabla["ha_mascara_agricola"] > 0, "agricola", "distrito")

censo = pd.read_parquet(DATOS / "cenagro" / "rec02_cultivos.parquet", columns=["ubigeo", "P025", "P026"])
censo["ubigeo"] = censo["ubigeo"].replace({"160109": "160801", "160114": "160803"})
censo = censo[censo["ubigeo"] != "120699"]
censo["ha_riego"] = censo["P025"].where(censo["P026"] == 1, 0.0)
riego = censo.groupby("ubigeo")[["ha_riego", "P025"]].sum()
tabla["fraccion_riego"] = (riego["ha_riego"] / riego["P025"]).reindex(tabla.index)

tabla = tabla.reset_index()
tabla.to_parquet(DATOS / "tierra_controles_distrito.parquet", index=False)
print(tabla.shape)
print(tabla.drop(columns=["ubigeo", "base_agregacion"]).describe().T[["count", "mean", "min", "max"]].round(2))
print("ha de la máscara agrícola:", round(tabla["ha_mascara_agricola"].sum()))
