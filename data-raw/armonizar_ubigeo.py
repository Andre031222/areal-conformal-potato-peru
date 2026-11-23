"""Asigna cada polígono distrital del INEI 2023 a un ubigeo del CENAGRO 2012.

Los ubigeos que existen en el censo se conservan (con las dos recodificaciones de Putumayo).
Cada distrito creado después de 2012 se asigna al distrito censal vecino de la misma provincia
con el que comparte más límite; si no tiene vecino censal en la provincia, al de cualquier
provincia. Se itera porque algunos distritos nuevos solo limitan con otros nuevos.
La asignación es una heurística: revisar `ubigeo_2023_a_2012.csv` contra las leyes de creación.

Salida: datos/ubigeo_2023_a_2012.csv y datos/limites/distritos_2012.gpkg (polígonos disueltos).
"""
import pathlib

import geopandas as gpd
import pandas as pd

DATOS = pathlib.Path(__file__).resolve().parents[1] / "data"
RECODIFICA = {"160109": "160801", "160114": "160803"}

censo = set(pd.read_csv(DATOS / "etiqueta_papa_distrito_2012.csv", dtype={"ubigeo": str})["ubigeo"])
pol = gpd.read_file(DATOS / "limites" / "DISTRITO.gpkg")
pol["ubigeo"] = pol["ubigeo"].astype(str).str.zfill(6)
pol = pol.dissolve(by="ubigeo", as_index=False, aggfunc="first")
pol["provincia"] = pol["ubigeo"].str[:4]

asignado = {u: u for u in pol["ubigeo"] if u in censo}
faltan_censo = censo - set(asignado)
print("polígonos:", len(pol), "| coinciden con el censo:", len(asignado), "| censo sin polígono:", sorted(faltan_censo))

metrico = pol.to_crs(32718).set_index("ubigeo")
borde = metrico.geometry.boundary
vecinos = gpd.sjoin(metrico[["geometry"]].reset_index(), metrico[["geometry"]].reset_index(),
                    predicate="touches")[["ubigeo_left", "ubigeo_right"]]
vecinos = vecinos[vecinos["ubigeo_left"] != vecinos["ubigeo_right"]]
vecinos["largo_m"] = [borde[a].intersection(borde[b]).length
                      for a, b in zip(vecinos["ubigeo_left"], vecinos["ubigeo_right"])]
vecinos = vecinos[vecinos["largo_m"] > 0]

filas = []
for ronda in range(1, 10):
    nuevos = [u for u in pol["ubigeo"] if u not in asignado]
    if not nuevos:
        break
    hechos = 0
    for u in nuevos:
        v = vecinos[(vecinos["ubigeo_left"] == u) & vecinos["ubigeo_right"].isin(asignado)]
        if v.empty:
            continue
        misma = v[v["ubigeo_right"].str[:4] == u[:4]]
        elegido = (misma if not misma.empty else v).sort_values("largo_m").iloc[-1]
        asignado[u] = asignado[elegido["ubigeo_right"]]
        filas.append({"ubigeo_2023": u, "ubigeo_2012": asignado[u], "via": elegido["ubigeo_right"],
                      "misma_provincia": not misma.empty, "borde_km": round(elegido["largo_m"] / 1000, 2),
                      "ronda": ronda})
        hechos += 1
    print("ronda", ronda, ":", hechos, "asignados")

sin = [u for u in pol["ubigeo"] if u not in asignado]
assert not sin, f"sin asignar: {sin}"
assert not faltan_censo, f"ubigeos del censo sin polígono: {faltan_censo}"

tabla = pd.DataFrame(filas)
nombres = pol.set_index("ubigeo")["nombdist"]
tabla["nombre_2023"] = tabla["ubigeo_2023"].map(nombres)
tabla["nombre_2012"] = tabla["ubigeo_2012"].map(nombres)
completo = pd.concat([
    pd.DataFrame({"ubigeo_2023": [u for u, d in asignado.items() if u == d],
                  "ubigeo_2012": [u for u, d in asignado.items() if u == d], "via": "", "ronda": 0}),
    tabla,
]).sort_values("ubigeo_2023")
completo.to_csv(DATOS / "ubigeo_2023_a_2012.csv", index=False)

pol["ubigeo_2012"] = pol["ubigeo"].map(asignado)
d12 = pol.dissolve(by="ubigeo_2012", as_index=False)[["ubigeo_2012", "geometry"]].rename(columns={"ubigeo_2012": "ubigeo"})
d12.to_file(DATOS / "limites" / "distritos_2012.gpkg", driver="GPKG")
print("distritos 2012:", len(d12), "| nuevos disueltos:", len(tabla),
      "| fuera de su provincia:", int((~tabla["misma_provincia"]).sum()))
print("recodificados de Putumayo presentes:", all(v in asignado for v in RECODIFICA.values()))
