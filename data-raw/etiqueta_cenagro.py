import pathlib

import pandas as pd

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
PAPA = [2610, 2611, 2612, 2613, 2614, 2615]
TAU = [0.005, 0.01, 0.05]

d = pd.read_parquet(DATOS / "cenagro" / "rec02_cultivos.parquet")
d["ubigeo"] = d["ubigeo"].replace({"160109": "160801", "160114": "160803"})
d = d[d["ubigeo"] != "120699"]
d["papa"] = d["P024_03"].isin(PAPA)
d["ha_papa"] = d["P025"].where(d["papa"], 0.0)
d["ha_papa_secano"] = d["P025"].where(d["papa"] & (d["P026"] == 2), 0.0)

assert abs(d["ha_papa"].sum() + 34.7 - 367_692) < 1, "los códigos de papa no reproducen el total oficial"

etiqueta = d.groupby("ubigeo").agg(
    ha_sembrada=("P025", "sum"),
    ha_papa=("ha_papa", "sum"),
    ha_papa_secano=("ha_papa_secano", "sum"),
    registros_cultivo=("papa", "size"),
    registros_papa=("papa", "sum"),
).reset_index()
etiqueta["fraccion_papa"] = etiqueta["ha_papa"] / etiqueta["ha_sembrada"]
etiqueta["fraccion_registros"] = etiqueta["registros_papa"] / etiqueta["registros_cultivo"]
for tau in TAU:
    etiqueta[f"y_tau{tau * 100:g}"] = (etiqueta["fraccion_papa"] >= tau).astype(int)

etiqueta.round(4).to_csv(DATOS / "etiqueta_papa_distrito_2012.csv", index=False)
print(len(etiqueta), "distritos;", "ha papa:", round(etiqueta["ha_papa"].sum(), 1))
print(etiqueta[[c for c in etiqueta if c.startswith("y_")]].mean().round(3).to_dict())
