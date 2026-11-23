#!/bin/sh
# Lanza las sensibilidades (20 pliegues por caso) y la validación con vecinos excluidos
# (1 783 distritos en tramos de 50), cada tarea en su propio proceso de R.
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p 04-resultados/sensibilidad/registros
tareas() {
  for caso in tau0.5 tau5 altitud_ph superficie pisco iid; do for f in $(seq 1 20); do echo "caso=$caso pliegue=$f"; done; done
  for i in $(seq 1 50 1783); do echo "caso=vecinos desde=$i hasta=$((i + 49))"; done
}
for p in "${PARALELO:-6}" 2 1; do
  tareas | xargs -P "$p" -L 1 sh -c 'Rscript 03-codigo/07_sensibilidad.R "$@" >> 04-resultados/sensibilidad/registros/$(echo "$@" | tr " =" "__").log 2>&1' _
  echo "pasada con $p procesos: $(find 04-resultados/sensibilidad -name "*.rds" | wc -l) de 156"
  [ "$(find 04-resultados/sensibilidad -name "*.rds" | wc -l)" -ge 156 ] && break
done
