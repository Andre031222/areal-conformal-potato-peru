"""Métricas fenológicas de NDVI por distrito 2012 (MODIS Terra y Aqua sobre la máscara agrícola).

1. Los compuestos de 16 días de los distritos 2023 se juntan en su distrito 2012, ponderando
   por la fracción agrícola (columna `peso`).
2. Serie mensual: máximo de los compuestos del mes; meses vacíos por interpolación lineal;
   suavizado Savitzky–Golay (ventana 5, grado 2).
3. Año fenológico propio de cada distrito: empieza en el mes de su NDVI mínimo climatológico
   (como en TIMESAT), para no cortar la estación verde; en la sierra suele empezar entre
   agosto y octubre. Se forman 24 años completos; cada uno se rotula con el año calendario
   de su séptimo mes. En cada año: NDVI máximo, mínimo, amplitud, mes del pico, inicio y fin
   de la estación verde (primer y último mes con NDVI ≥ mínimo + 20 % de la amplitud) y su
   duración en meses. Los meses del pico, inicio y fin son calendario (1 = enero).

Salida: datos/fenologia_anual_distrito.parquet (ubigeo × año agrícola).
"""
import pathlib

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
UMBRAL = 0.20

ndvi = pd.read_parquet(DATOS / "modis" / "ndvi_16dias_distrito.parquet")
mapa = pd.read_csv(DATOS / "ubigeo_2023_a_2012.csv", dtype=str).set_index("ubigeo_2023")["ubigeo_2012"]
ndvi["ubigeo"] = ndvi["ubigeo"].astype(str).str.zfill(6)
sin_mapa = set(ndvi["ubigeo"]) - set(mapa.index)
print("ubigeos NDVI sin correspondencia 2012 (se descartan):", sorted(sin_mapa))
ndvi = ndvi[ndvi["ubigeo"].isin(mapa.index)]
ndvi["ubigeo"] = ndvi["ubigeo"].map(mapa)
ndvi["suma"] = ndvi["ndvi"] * ndvi["peso"]
compuesto = ndvi.groupby(["ubigeo", "fecha"])[["suma", "peso"]].sum()
compuesto = (compuesto["suma"] / compuesto["peso"]).rename("ndvi").reset_index()

compuesto["mes"] = compuesto["fecha"].dt.to_period("M")
mensual = compuesto.groupby(["ubigeo", "mes"])["ndvi"].max().unstack("mes")
meses = pd.period_range("2000-07", "2025-06", freq="M")
mensual = mensual.reindex(columns=meses)
vacios = mensual.isna().mean(axis=1)
mensual = mensual.T.interpolate(limit_direction="both").T
suave = savgol_filter(mensual.to_numpy(), window_length=5, polyorder=2, axis=1)

climatologia = np.stack([suave[:, m::12].mean(1) for m in range(12)], axis=1)
desfase = climatologia.argmin(1)
mes_inicio = (desfase + 6) % 12 + 1
filas = []
for k in range(24):
    j = desfase[:, None] + 12 * k + np.arange(12)[None, :]
    x = np.take_along_axis(suave, j, axis=1)
    vmax, vmin = x.max(1), x.min(1)
    amp = vmax - vmin
    verde = x >= (vmin + UMBRAL * amp)[:, None]
    inicio = verde.argmax(1)
    fin = 11 - verde[:, ::-1].argmax(1)
    calendario = lambda i: (mes_inicio - 1 + i) % 12 + 1
    filas.append(pd.DataFrame({
        "ubigeo": mensual.index, "anio_agricola": meses[j[:, 6]].year,
        "mes_inicio_anio": mes_inicio,
        "ndvi_max": vmax, "ndvi_min": vmin, "ndvi_amplitud": amp,
        "mes_pico": calendario(x.argmax(1)), "inicio_verde": calendario(inicio),
        "fin_verde": calendario(fin), "duracion_verde_meses": fin - inicio + 1,
        "fraccion_meses_vacios": vacios.to_numpy(),
    }))

panel = pd.concat(filas, ignore_index=True)
panel.to_parquet(DATOS / "fenologia_anual_distrito.parquet", index=False)
print(panel.shape, "| distritos:", panel["ubigeo"].nunique())
print(panel.drop(columns=["ubigeo", "anio_agricola"]).describe().T[["mean", "min", "max"]].round(3))
