#!/bin/sh
# Construye data/tabla_distrital_2012.csv desde las fuentes públicas ya descargadas en data/
# (descarga: sh data-raw/descargar.sh todo; ver data-raw/FUENTES.md). Correr desde la raíz del repositorio.
set -e
cd "$(dirname "$0")/.."
for paso in etiqueta_cenagro armonizar_ubigeo puntos_agricolas clima_terraclimate tierra_controles fenologia_ndvi tabla_distrital; do
  echo "== $paso"
  python3 -W ignore "data-raw/$paso.py"
done
