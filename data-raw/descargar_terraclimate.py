import pathlib
import socket
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import xarray as xr

VARIABLES = ["tmax", "tmin", "ppt", "pet", "def", "soil"]
ANIOS = range(1981, 2026)
URL = ("http://thredds.northwestknowledge.net:8080/thredds/ncss/TERRACLIMATE_ALL/data/"
       "TerraClimate_{v}_{a}.nc?var={v}&north=0.5&south=-18.6&west=-81.5&east=-68.5"
       "&temporal=all&accept=netcdf")
RAIZ = pathlib.Path(__file__).resolve().parents[1] / "data" / "terraclimate"
socket.setdefaulttimeout(180)


def completo(ruta: pathlib.Path) -> bool:
    if not ruta.exists():
        return False
    try:
        with xr.open_dataset(ruta) as ds:
            return ds.sizes.get("time", 0) == 12
    except Exception:
        return False


def bajar(v: str, a: int) -> str:
    destino = RAIZ / v / f"TerraClimate_{v}_{a}_peru.nc"
    if completo(destino):
        return f"ya estaba {v} {a}"
    temporal = destino.with_suffix(".parcial")
    for intento in range(5):
        try:
            urllib.request.urlretrieve(URL.format(v=v, a=a), temporal)
            temporal.rename(destino)
            if completo(destino):
                return f"listo {v} {a}"
        except Exception as error:
            print("reintento", v, a, error, flush=True)
        time.sleep(10 * (intento + 1))
    return f"FALLÓ {v} {a}"


for v in VARIABLES:
    (RAIZ / v).mkdir(parents=True, exist_ok=True)
with ThreadPoolExecutor(max_workers=4) as grupo:
    for mensaje in grupo.map(lambda p: bajar(*p), [(v, a) for v in VARIABLES for a in ANIOS]):
        print(mensaje, flush=True)
