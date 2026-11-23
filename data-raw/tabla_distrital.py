"""Tabla distrital final (etapa 2): una fila por distrito del IV CENAGRO 2012.

Periodos del clima (campaña c = septiembre de c-1 a mayo de c):
- `_ajuste`: campañas 1991–2012, las dos décadas previas al censo (ajuste del modelo).
- `_p1`: campañas 1982–2000 (referencia; 1981 no tiene septiembre-diciembre de 1980).
- `_p2`: campañas 2001–2025 (reciente).
Fenología: promedio de los años 2001–2012; el mes del pico es la mediana.
Sierra: altitud media de la superficie agrícola ≥ 2 300 m.

Salida: datos/tabla_distrital_2012.csv y .parquet. Corre las validaciones antes de guardar.
"""
import pathlib

import numpy as np
import pandas as pd

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
PERIODOS = {"ajuste": (1991, 2012), "p1": (1982, 2000), "p2": (2001, 2025)}
SIERRA_M = 2300

etiqueta = pd.read_csv(DATOS / "etiqueta_papa_distrito_2012.csv", dtype={"ubigeo": str})
tierra = pd.read_parquet(DATOS / "tierra_controles_distrito.parquet")
clima = pd.read_parquet(DATOS / "clima_campana_distrito.parquet")
feno = pd.read_parquet(DATOS / "fenologia_anual_distrito.parquet")

variables_clima = [c for c in clima.columns if c not in ("ubigeo", "campana")]
tabla = etiqueta.merge(tierra, on="ubigeo", how="left", validate="one_to_one")
for nombre, (a, b) in PERIODOS.items():
    p = clima[clima["campana"].between(a, b)].groupby("ubigeo")[variables_clima].mean()
    assert clima[clima["campana"].between(a, b)]["campana"].nunique() == b - a + 1
    tabla = tabla.merge(p.add_suffix(f"_{nombre}"), left_on="ubigeo", right_index=True, how="left")

f = feno[feno["anio_agricola"].between(2001, 2012)].groupby("ubigeo")
fen = f[["ndvi_max", "ndvi_min", "ndvi_amplitud", "duracion_verde_meses"]].mean()
fen["mes_pico"] = f["mes_pico"].median()
fen["inicio_verde"] = f["inicio_verde"].median()
fen["mes_inicio_anio"] = f["mes_inicio_anio"].first()
tabla = tabla.merge(fen.add_prefix("feno_"), left_on="ubigeo", right_index=True, how="left")

tabla["sierra"] = (tabla["altitud_media_m"] >= SIERRA_M).astype(int)
tabla.insert(1, "departamento", tabla["ubigeo"].str[:2])
tabla.insert(2, "provincia", tabla["ubigeo"].str[:4])

# Validaciones
assert len(tabla) == 1783 and tabla["ubigeo"].is_unique
assert abs(tabla["ha_papa"].sum() + 34.7 - 367_692) < 1, "la papa no suma el total oficial"
assert tabla[[f"ppt_mm_{p}" for p in PERIODOS]].notna().all().all(), "falta clima en algún distrito"
assert tabla["altitud_media_m"].between(-10, 5500).all()
assert tabla["tmin_c_ajuste"].lt(tabla["tmax_c_ajuste"]).all()
assert tabla["ph"].dropna().between(3.5, 9.5).all()
r_alt_tmin = tabla["altitud_media_m"].corr(tabla["tmin_c_ajuste"])
assert r_alt_tmin < -0.8, f"altitud y Tmín deberían correlacionar fuerte en negativo ({r_alt_tmin:.2f})"

faltantes = tabla.isna().sum()
print("distritos:", len(tabla), "| columnas:", tabla.shape[1])
print("faltantes por columna (solo las que tienen):")
print(faltantes[faltantes > 0].to_string())
print("sin suelo:", tabla.loc[tabla["ph"].isna(), "ubigeo"].tolist())
print("sin fenología:", int(tabla["feno_ndvi_max"].isna().sum()), "distritos")
print("sierra:", int(tabla["sierra"].sum()), "| con papa (τ = 1 %) en la sierra:",
      int(tabla.loc[tabla["sierra"] == 1, "y_tau1"].sum()))
print("corr(altitud, Tmín ajuste):", round(r_alt_tmin, 3))
cambio = (tabla["tmin_c_p2"] - tabla["tmin_c_p1"]).describe().round(3)
print("cambio de Tmín de campaña p2 - p1 (°C):", cambio[["mean", "min", "max"]].to_dict())
print("cambio de meses con helada en la sierra p2 - p1:",
      round((tabla["meses_helada_p2"] - tabla["meses_helada_p1"])[tabla["sierra"] == 1].mean(), 3))

tabla.to_parquet(DATOS / "tabla_distrital_2012.parquet", index=False)
tabla.round(4).to_csv(DATOS / "tabla_distrital_2012.csv", index=False)
print("guardado: datos/tabla_distrital_2012.csv y .parquet")
