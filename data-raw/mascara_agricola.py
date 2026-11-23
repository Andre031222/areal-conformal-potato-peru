import pathlib

import numpy as np
import pyogrio
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
FUENTE = DATOS / "superficie_agricola" / "MNSA_2024" / "SuperficieAgricolaNacional.shp"
SALIDA = DATOS / "superficie_agricola" / "mascara_agricola_2024.tif"
PASO = 0.00125
OESTE, NORTE, ESTE, SUR = -81.25, 0.0, -68.75, -18.5

ancho = round((ESTE - OESTE) / PASO)
alto = round((NORTE - SUR) / PASO)
transformada = from_origin(OESTE, NORTE, PASO, PASO)
mascara = np.zeros((alto, ancho), dtype="uint8")

total = pyogrio.read_info(FUENTE)["features"]
for inicio in range(0, total, 50):
    lote = pyogrio.read_dataframe(FUENTE, columns=[], skip_features=inicio, max_features=50)
    rasterize(((g, 1) for g in lote.geometry if g is not None), out=mascara, transform=transformada)
    print(inicio + len(lote), "de", total, flush=True)

perfil = {"driver": "GTiff", "height": alto, "width": ancho, "count": 1, "dtype": "uint8",
          "crs": "EPSG:4326", "transform": transformada, "compress": "deflate", "tiled": True, "nodata": None}
with rasterio.open(SALIDA, "w", **perfil) as destino:
    destino.write(mascara, 1)
print("hectáreas aproximadas:", round(mascara.sum() * (PASO * 111_320) ** 2 / 10_000))
