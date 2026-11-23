"""Píxeles de la máscara agrícola MIDAGRI 2024 (0,00125°) con su distrito 2012.

Cada píxel agrícola es un punto donde se leen todas las capas (clima, suelo, relieve, acceso).
Los distritos sin píxeles agrícolas usan todos sus píxeles con un paso de 4, marcados con
base = "distrito", para no perderlos de la tabla.

Salida: datos/puntos_agricolas.parquet (lon, lat, ubigeo, base).
"""
import pathlib

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
PASO_RESPALDO = 4

d12 = gpd.read_file(DATOS / "limites" / "distritos_2012.gpkg")
d12["id"] = np.arange(1, len(d12) + 1, dtype=np.int16)

with rasterio.open(DATOS / "superficie_agricola" / "mascara_agricola_2024.tif") as r:
    mascara = r.read(1).astype(bool)
    transformada = r.transform

ids = rasterize(zip(d12.geometry, d12["id"]), out_shape=mascara.shape, transform=transformada,
                fill=0, dtype="int16")

fila, col = np.nonzero(mascara & (ids > 0))
base = np.zeros(len(fila), dtype=bool)
con_agro = np.unique(ids[fila, col])
sin_agro = np.setdiff1d(d12["id"].to_numpy(), con_agro)
if len(sin_agro):
    sub = ids[::PASO_RESPALDO, ::PASO_RESPALDO]
    f2, c2 = np.nonzero(np.isin(sub, sin_agro))
    f2, c2 = f2 * PASO_RESPALDO, c2 * PASO_RESPALDO
    fila, col = np.concatenate([fila, f2]), np.concatenate([col, c2])
    base = np.concatenate([base, np.ones(len(f2), dtype=bool)])

lon = transformada.c + (col + 0.5) * transformada.a
lat = transformada.f + (fila + 0.5) * transformada.e
ubigeo = d12["ubigeo"].to_numpy()[ids[fila, col] - 1]

puntos = pd.DataFrame({"lon": lon, "lat": lat, "ubigeo": ubigeo,
                       "base": np.where(base, "distrito", "agricola")})
puntos.to_parquet(DATOS / "puntos_agricolas.parquet", index=False)

fuera = int((mascara & (ids == 0)).sum())
print("píxeles agrícolas:", int((~base).sum()), "| fuera de todo distrito:", fuera)
print("distritos con píxeles agrícolas:", len(con_agro), "de", len(d12),
      "| sin píxeles agrícolas (respaldo con todo el distrito):", len(sin_agro))
faltan = set(d12["ubigeo"]) - set(puntos["ubigeo"])
print("distritos sin ningún punto:", sorted(faltan))
