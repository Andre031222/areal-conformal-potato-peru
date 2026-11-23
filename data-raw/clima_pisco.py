"""Clima de la campaña agrícola por distrito con PISCO (sensibilidad 7.4 del plan de análisis).

Mismas métricas y mismo método que clima_terraclimate.py, con las grillas PISCO del SENAMHI
(Aybar et al., 2020; Huerta et al., 2023), construidas con la red nacional de estaciones:
PISCOp v2.1 (precipitación) y PISCOt v1.1 (Tmáx y Tmín), mensuales, versión estable, 0,1°,
enero de 1981 a diciembre de 2016.

- Las métricas se calculan en cada celda y se promedian en el distrito ponderando por el número
  de píxeles agrícolas de la celda (datos/puntos_agricolas.parquet).
- Campaña c = septiembre de c-1 a mayo de c; periodo de ajuste = campañas 1991-2012, como en la
  tabla principal.
- PISCO no trae evapotranspiración: ET0 por Hargreaves-Samani (FAO-56, ec. 52) con Tmáx, Tmín y la
  radiación extraterrestre del día 15 de cada mes (FAO-56, ec. 21).
- SPEI-3 como en clima_terraclimate.py (logística generalizada por L-momentos), con 1981-2016
  como referencia.

Salida: datos/clima_pisco_distrito_2012.csv (ubigeo y las cuatro variables climáticas de M1).
"""
import pathlib

import numpy as np
import pandas as pd
import xarray as xr
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import norm

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
PISCO = DATOS / "pisco"
MESES = 12 * (2016 - 1981 + 1)


def leer(archivo, variable):
    x = xr.open_dataset(PISCO / archivo, decode_times=False)[variable]
    t = x["T"].to_numpy()
    assert x.sizes["T"] == MESES and t[0] == 252.5, f"{archivo}: se esperaba ene-1981 a dic-2016"
    return x.transpose("T", "Y", "X")


def celdas_de(x, puntos):
    """Celda de la grilla de cada punto; los puntos en celdas sin dato van a la celda válida más cercana."""
    lat, lon = x["Y"].to_numpy(), x["X"].to_numpy()
    fila = np.clip(np.rint((puntos["lat"].to_numpy() - lat[0]) / (lat[1] - lat[0])).astype(int), 0, len(lat) - 1)
    col = np.clip(np.rint((puntos["lon"].to_numpy() - lon[0]) / (lon[1] - lon[0])).astype(int), 0, len(lon) - 1)
    celda = fila * len(lon) + col
    valida = np.isfinite(x.isel(T=slice(0, 24)).to_numpy()).all(axis=0).ravel()
    malas = ~valida[celda]
    if malas.any():
        fv, cv = np.divmod(np.nonzero(valida)[0], len(lon))
        _, cerca = cKDTree(np.column_stack([fv, cv])).query(np.column_stack([fila[malas], col[malas]]))
        celda[malas] = np.nonzero(valida)[0][cerca]
    print(f"  puntos reasignados a la celda válida más cercana: {int(malas.sum())}")
    return celda, lat


def spei3(balance):
    """Igual que en clima_terraclimate.py. Logística generalizada por L-momentos insesgados (Hosking, 1997), como el paquete SPEI de R."""
    acumulado = pd.DataFrame(balance).rolling(3).sum().to_numpy()
    salida = np.full_like(acumulado, np.nan)
    for mes in range(12):
        filas = np.arange(mes, len(acumulado), 12)
        filas = filas[np.isfinite(acumulado[filas, 0])]
        x = np.sort(acumulado[filas], axis=0)
        n = len(x)
        i = np.arange(n)[:, None]
        b0 = x.mean(0)
        b1 = (i / (n - 1) * x).mean(0)
        b2 = (i * (i - 1) / ((n - 1) * (n - 2)) * x).mean(0)
        l1, l2, l3 = b0, 2 * b1 - b0, 6 * b2 - 6 * b1 + b0
        k = -l3 / l2
        k = np.where(np.abs(k) < 1e-6, 1e-6, k)
        alfa = l2 * np.sin(k * np.pi) / (k * np.pi)
        xi = l1 - alfa * (1 / k - np.pi / np.sin(k * np.pi))
        arg = 1 - k * (acumulado[filas] - xi) / alfa
        with np.errstate(divide="ignore", invalid="ignore"):
            y = -np.log(arg) / k
        y = np.where(arg > 0, y, np.where(k > 0, np.inf, -np.inf))
        prob = 1 / (1 + np.exp(-y))
        salida[filas] = norm.ppf(np.clip(prob, 1e-6, 1 - 1e-6))
    return salida


def ra_mm_dia(lat_grados, mes):
    """Radiación extraterrestre (mm/día equivalentes) el día 15 del mes, FAO-56 ec. 21."""
    j = pd.Timestamp(2001, mes, 15).dayofyear
    phi = np.deg2rad(lat_grados)
    dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365)
    delta = 0.409 * np.sin(2 * np.pi * j / 365 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1, 1))
    ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(ws))
    return 0.408 * ra


puntos = pd.read_parquet(DATOS / "puntos_agricolas.parquet")
distritos, distrito_local = np.unique(puntos["ubigeo"].to_numpy(), return_inverse=True)
fechas = pd.date_range("1981-01-01", periods=MESES, freq="MS")
dias = fechas.days_in_month.to_numpy()[:, None]

prec = leer("PISCOp_v2.1_stable_mensual.nc", "Prec")
tmin = leer("PISCOt_v1.1_tmin_stable_mensual.nc", "tmin")
tmax = leer("PISCOt_v1.1_tmax_stable_mensual.nc", "tmax")
assert np.allclose(tmin["Y"], tmax["Y"]) and np.allclose(tmin["X"], tmax["X"])


def matriz(x, celda):
    usadas, local = np.unique(celda, return_inverse=True)
    W = csr_matrix((np.ones(len(local)), (distrito_local, local)), shape=(len(distritos), len(usadas)))
    W = csr_matrix(W.multiply(1.0 / np.asarray(W.sum(axis=1))))
    valores = x.to_numpy().reshape(MESES, -1)[:, usadas].astype("float64")
    return W, usadas, valores


print("precipitación:")
celda_p, _ = celdas_de(prec, puntos)
Wp, usadas_p, ppt = matriz(prec, celda_p)
print("temperatura:")
celda_t, lat_t = celdas_de(tmin, puntos)
Wt, usadas_t, tn = matriz(tmin, celda_t)
_, _, tx = matriz(tmax, celda_t)
for nombre, v in {"ppt": ppt, "tmin": tn, "tmax": tx}.items():
    assert np.isfinite(v).all(), f"{nombre}: faltan valores en celdas usadas"
print("celdas usadas: precipitación", len(usadas_p), "| temperatura", len(usadas_t))

# ET0 de Hargreaves en las celdas de temperatura; el balance P - ET0 necesita ambas en la misma
# celda, así que se calcula por punto agrícola único (par de celdas) y se agrega igual que lo demás.
lat_celda = lat_t[np.divmod(usadas_t, tmin.sizes["X"])[0]]
ra = np.stack([ra_mm_dia(lat_celda, m) for m in fechas.month])
et0 = 0.0023 * ra * ((tx + tn) / 2 + 17.8) * np.sqrt(np.maximum(tx - tn, 0)) * dias

pares, par_local = np.unique(np.column_stack([celda_p, celda_t]), axis=0, return_inverse=True)
par_local = par_local.ravel()
ip = np.searchsorted(usadas_p, pares[:, 0]); it = np.searchsorted(usadas_t, pares[:, 1])
Wb = csr_matrix((np.ones(len(par_local)), (distrito_local, par_local)), shape=(len(distritos), len(pares)))
Wb = csr_matrix(Wb.multiply(1.0 / np.asarray(Wb.sum(axis=1))))
print("pares de celdas (precipitación, temperatura):", len(pares))

spei = spei3(ppt[:, ip] - et0[:, it])

filas = []
for c in range(1991, 2013):
    t = np.nonzero((fechas >= f"{c - 1}-09-01") & (fechas <= f"{c}-05-01"))[0]
    assert len(t) == 9
    filas.append(pd.DataFrame({
        "ubigeo": distritos, "campana": c,
        "ppt_mm": Wp @ ppt[t].sum(0),
        "tmin_c": Wt @ tn[t].mean(0),
        "meses_helada": Wt @ (tn[t] <= 0).sum(0),
        "spei3_min": Wb @ spei[t].min(0),
    }))
panel = pd.concat(filas, ignore_index=True)
tabla = panel.groupby("ubigeo")[["ppt_mm", "tmin_c", "meses_helada", "spei3_min"]].mean().add_suffix("_pisco").reset_index()
assert len(tabla) == 1783 and tabla.notna().all().all()
tabla.round(4).to_csv(DATOS / "clima_pisco_distrito_2012.csv", index=False)

tc = pd.read_csv(DATOS / "tabla_distrital_2012.csv", dtype={"ubigeo": str}).set_index("ubigeo")
tabla = tabla.set_index("ubigeo")
print("\nPISCO frente a TerraClimate, campañas 1991-2012 (n = 1783):")
for v, w in [("tmin_c", "tmin_c_ajuste"), ("ppt_mm", "ppt_mm_ajuste"), ("meses_helada", "meses_helada_ajuste"), ("spei3_min", "spei3_min_ajuste")]:
    a, b = tabla[f"{v}_pisco"], tc.loc[tabla.index, w]
    print(f"  {v:13s} r = {a.corr(b):.3f} | media PISCO {a.mean():8.2f} | media TerraClimate {b.mean():8.2f}")
