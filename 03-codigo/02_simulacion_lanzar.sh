#!/bin/sh
# Lanza las 60 tareas de la simulación (3 valores de phi × 20 réplicas), PARALELO a la vez (8 por defecto), cada una
# en un proceso de R independiente. Las tareas ya guardadas en partes_v2 se saltan.
# Correr desde la raíz del repositorio: sh 03-codigo/02_simulacion_lanzar.sh
cd "$(dirname "$0")/.."
# Un solo hilo de cálculo por proceso: el paralelismo lo da xargs (si no, BLAS usa todos los núcleos
# en cada proceso, la carga se dispara y el equipo se calienta).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p 04-resultados/simulacion/registros
tareas() { for phi in 0.9 0.5 0.2; do for r in $(seq 1 20); do echo "$phi $r"; done; done; }
# INLA falla de forma intermitente cuando la máquina está muy cargada. Primera pasada con
# PARALELO procesos; luego dos pasadas más con menos carga completan lo que falte (las tareas
# ya guardadas se saltan).
for p in "${PARALELO:-8}" 4 1; do
  tareas | xargs -P "$p" -n 2 sh -c 'Rscript 03-codigo/02_simulacion.R phi=$0 rep=$1 >> 04-resultados/simulacion/registros/phi$0_rep$1.log 2>&1'
  echo "pasada con $p procesos: $(ls 04-resultados/simulacion/partes_v2 | wc -l) de 60"
  [ "$(ls 04-resultados/simulacion/partes_v2 | wc -l)" -ge 60 ] && break
done
echo "simulación terminada: $(ls 04-resultados/simulacion/partes_v2 | wc -l) de 60 tareas"
