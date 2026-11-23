#!/bin/sh
# Corre los 20 pliegues de la aplicación al censo, cada uno en su propio proceso de R.
# Los pliegues ya guardados se saltan; las pasadas con menos procesos rellenan los que fallen.
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p 04-resultados/aplicacion/registros 04-resultados/aplicacion/pliegues
for p in "${PARALELO:-4}" 2 1; do
  seq 1 20 | xargs -P "$p" -I{} sh -c 'Rscript 03-codigo/04_aplicacion_censo.R pliegue={} >> 04-resultados/aplicacion/registros/pliegue{}.log 2>&1'
  echo "pasada con $p procesos: $(ls 04-resultados/aplicacion/pliegues | wc -l) de 20"
  [ "$(ls 04-resultados/aplicacion/pliegues | wc -l)" -ge 20 ] && break
done
