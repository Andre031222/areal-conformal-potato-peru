import pathlib
import urllib.request
import zipfile

import pandas as pd
import pyreadstat

RAIZ = pathlib.Path(__file__).resolve().parents[1] / "data" / "cenagro" / "spss"
RAIZ.mkdir(parents=True, exist_ok=True)
URL = "https://proyectos.inei.gob.pe/iinei/srienaho/descarga/SPSS/{c}-Modulo231.zip"
COLUMNAS = ["P001", "P002", "P003", "NPRIN", "NPARCX", "P024_01", "P024_03", "P025", "P026"]

tablas = []
for codigo in range(337, 362):
    destino = RAIZ / f"{codigo}-Modulo231.zip"
    if not destino.exists() or destino.stat().st_size < 100_000:
        urllib.request.urlretrieve(URL.format(c=codigo), destino)
    with zipfile.ZipFile(destino) as paquete:
        nombre = next(n for n in paquete.namelist() if n.upper().endswith("REC02.SAV"))
        paquete.extract(nombre, RAIZ / str(codigo))
    df, _ = pyreadstat.read_sav(str(RAIZ / str(codigo) / nombre), usecols=COLUMNAS)
    df = df[df["P024_01"] < 90]
    tablas.append(df)
    print(codigo, len(df), flush=True)

todo = pd.concat(tablas, ignore_index=True)
todo["ubigeo"] = todo["P001"].str.zfill(2) + todo["P002"].str.zfill(2) + todo["P003"].str.zfill(2)
todo.to_parquet(RAIZ.parent / "rec02_cultivos.parquet", index=False)
print("filas", len(todo))
