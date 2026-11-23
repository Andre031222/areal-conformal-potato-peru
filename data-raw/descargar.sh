#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../data"
QUE="${1:-limites}"

limites() {
  mkdir -p limites
  curl -fL --retry 3 -o limites/Distrito.rar https://ide.inei.gob.pe/files/Distrito.rar
}

cenagro() {
  mkdir -p cenagro
  curl -fL --retry 3 -C - -o cenagro/2012.zip \
    https://www.inei.gob.pe/media/DATOS_ABIERTOS/CENAGRO/DATA/2012.zip
}

enaho() {
  mkdir -p enaho
  for c in 498 546 603 634 687 737 759 784 906 966; do
    [ -s "enaho/$c-Modulo22.zip" ] || curl -fL --retry 3 -o "enaho/$c-Modulo22.zip" \
      "https://proyectos.inei.gob.pe/iinei/srienaho/descarga/STATA/$c-Modulo22.zip"
  done
}

terraclimate() {
  python3 ../data-raw/descargar_terraclimate.py
}

pisco() {
  mkdir -p pisco
  local base=https://iridl.ldeo.columbia.edu/SOURCES/.SENAMHI/.HSR/.PISCO/.Prec/.v2p1
  curl -fL -o pisco/PISCOp_v2.1_stable_mensual.nc "$base/.stable/.monthly/.Prec/data.nc"
  curl -fL -o pisco/PISCOp_v2.1_unstable_mensual.nc "$base/.unstable/.monthly/.Prec/data.nc"
  local temp=https://iridl.ldeo.columbia.edu/SOURCES/.SENAMHI/.HSR/.PISCO/.Temp/.v1p1
  for v in tmax tmin; do
    curl -fL -o "pisco/PISCOt_v1.1_${v}_stable_mensual.nc" "$temp/.$v/.stable/.monthly/.$v/data.nc"
  done
}

cenagro_spss() {
  python3 ../data-raw/descargar_cenagro_spss.py
}

suelo_relieve() {
  python3 ../data-raw/descargar_suelo_relieve.py
}

case "$QUE" in
  limites) limites ;;
  cenagro) cenagro ;;
  cenagro_spss) cenagro_spss ;;
  enaho) enaho ;;
  terraclimate) terraclimate ;;
  pisco) pisco ;;
  suelo_relieve) suelo_relieve ;;
  todo) limites; cenagro_spss; enaho; terraclimate; pisco; suelo_relieve ;;
  *) echo "Opción desconocida: $QUE" >&2; exit 1 ;;
esac
