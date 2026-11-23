"""Clima de la campaña agrícola por distrito 2012, TerraClimate 1981–2025.

Campaña c = septiembre del año c-1 a mayo del año c (ciclo de la papa altoandina).
Con datos desde enero de 1981, la primera campaña completa es 1982 y la última 2025.

Las métricas se calculan en cada celda de TerraClimate (1/24°) y después se promedian en el
distrito ponderando por el número de píxeles agrícolas de cada celda. Así las métricas no
lineales (heladas, grados-día, SPEI) se calculan donde ocurren y no sobre un promedio.

SPEI-3 (Vicente-Serrano et al., 2010): balance P - ET0 acumulado en 3 meses, distribución
logística generalizada ajustada por L-momentos insesgados, como el paquete SPEI de R
(Beguería et al., 2014), por celda y por mes calendario, con todo 1981–2025 como referencia.

Salida: datos/clima_campana_distrito.parquet (ubigeo × campaña).
"""
import pathlib

import numpy as np
import pandas as pd
import xarray as xr
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import norm

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
TC = DATOS / "terraclimate"
ANIOS = range(1981, 2026)
VARIABLES = ["ppt", "tmax", "tmin", "pet", "def", "soil"]
BASE_GDD = 5.0


def leer(variable):
    partes = [xr.open_dataset(TC / variable / f"TerraClimate_{variable}_{a}_peru.nc")[variable] for a in ANIOS]
    serie = xr.concat(partes, dim="time")
    assert serie.sizes["time"] == 12 * len(ANIOS), f"{variable}: faltan meses"
    return serie


puntos = pd.read_parquet(DATOS / "puntos_agricolas.parquet")
ppt = leer("ppt")
lat, lon = ppt["lat"].to_numpy(), ppt["lon"].to_numpy()
nlat, nlon = len(lat), len(lon)

fila = np.clip(np.rint((puntos["lat"].to_numpy() - lat[0]) / (lat[1] - lat[0])).astype(int), 0, nlat - 1)
col = np.clip(np.rint((puntos["lon"].to_numpy() - lon[0]) / (lon[1] - lon[0])).astype(int), 0, nlon - 1)
celda = fila * nlon + col

valida = np.isfinite(ppt.isel(time=0).to_numpy()).ravel()
malas = ~valida[celda]
if malas.any():
    fv, cv = np.divmod(np.nonzero(valida)[0], nlon)
    arbol = cKDTree(np.column_stack([fv, cv]))
    _, cerca = arbol.query(np.column_stack([fila[malas], col[malas]]))
    celda[malas] = np.nonzero(valida)[0][cerca]
print("puntos reasignados a la celda terrestre más cercana:", int(malas.sum()))

usadas, celda_local = np.unique(celda, return_inverse=True)
distritos, distrito_local = np.unique(puntos["ubigeo"].to_numpy(), return_inverse=True)
W = csr_matrix((np.ones(len(puntos)), (distrito_local, celda_local)), shape=(len(distritos), len(usadas)))
W = csr_matrix(W.multiply(1.0 / np.asarray(W.sum(axis=1))))
print("celdas TerraClimate usadas:", len(usadas), "| distritos:", len(distritos))

datos = {"ppt": ppt.to_numpy().reshape(len(ppt["time"]), -1)[:, usadas].astype("float64")}
for v in VARIABLES[1:]:
    datos[v] = leer(v).to_numpy().reshape(12 * len(ANIOS), -1)[:, usadas].astype("float64")
for v, x in datos.items():
    assert np.isfinite(x).all(), f"{v}: hay valores faltantes en celdas usadas"

fechas = pd.date_range("1981-01-01", periods=12 * len(ANIOS), freq="MS")
dias = fechas.days_in_month.to_numpy()[:, None]


def spei3(balance):
    """Logística generalizada por L-momentos insesgados (Hosking, 1997), como el paquete SPEI de R."""
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


def promedio(valores):
    """Promedio ponderado del distrito que ignora las celdas sin valor."""
    ok = np.isfinite(valores)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (W @ np.where(ok, valores, 0.0)) / (W @ ok.astype(float))


spei = spei3(datos["ppt"] - datos["pet"])
print("celdas-mes sin SPEI-3 (desde marzo de 1981):", int(np.isnan(spei[2:]).sum()))
gdd = np.maximum((datos["tmax"] + datos["tmin"]) / 2 - BASE_GDD, 0) * dias

filas = []
for c in range(1982, 2026):
    t = np.nonzero((fechas >= f"{c - 1}-09-01") & (fechas <= f"{c}-05-01"))[0]
    assert len(t) == 9
    celdas = {
        "ppt_mm": datos["ppt"][t].sum(0),
        "tmax_c": datos["tmax"][t].mean(0),
        "tmin_c": datos["tmin"][t].mean(0),
        "meses_helada": (datos["tmin"][t] <= 0).sum(0),
        "gdd5": gdd[t].sum(0),
        "et0_mm": datos["pet"][t].sum(0),
        "deficit_mm": datos["def"][t].sum(0),
        "humedad_suelo_mm": datos["soil"][t].mean(0),
        "spei3_medio": spei[t].mean(0),
        "spei3_min": spei[t].min(0),
    }
    fila_c = pd.DataFrame({k: promedio(v) for k, v in celdas.items()})
    fila_c.insert(0, "campana", c)
    fila_c.insert(0, "ubigeo", distritos)
    filas.append(fila_c)

panel = pd.concat(filas, ignore_index=True)
panel.to_parquet(DATOS / "clima_campana_distrito.parquet", index=False)
print(panel.shape, "| campañas", panel["campana"].min(), "-", panel["campana"].max())
print(panel.drop(columns=["ubigeo", "campana"]).describe().T[["mean", "min", "max"]].round(2))
