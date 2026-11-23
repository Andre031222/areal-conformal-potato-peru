#!/bin/sh
# Desviación 11: simulación con el modelo mal especificado (phi = 0,5; 20 réplicas por escenario),
# cada tarea en su propio proceso de R. Pasadas con menos procesos rellenan las que fallen.
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p 04-resultados/simulacion/registros
ESCENARIOS="${ESCENARIOS:-colas no_estacionario no_estacionario_este}"
TOTAL=$((20 * $(echo $ESCENARIOS | wc -w)))
tareas() { for e in $ESCENARIOS; do for r in $(seq 1 20); do echo "$e $r"; done; done; }
for p in "${PARALELO:-8}" 4 1; do
  tareas | xargs -P "$p" -n 2 sh -c 'Rscript 03-codigo/02_simulacion.R escenario=$0 phi=0.5 rep=$1 >> 04-resultados/simulacion/registros/mal_$0_rep$1.log 2>&1'
  echo "pasada con $p procesos: $(ls 04-resultados/simulacion/partes_mal | grep -cE "^($(echo $ESCENARIOS | tr ' ' '|'))_phi") de $TOTAL"
  [ "$(ls 04-resultados/simulacion/partes_mal | grep -cE "^($(echo $ESCENARIOS | tr ' ' '|'))_phi")" -ge $TOTAL ] && break
done
