#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=4
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
mkdir -p hybrid
sha256sum src/SingleIndexGraph/*.jl experiments/single_index/mave_hybrid* experiments/single_index/mave_confirmation.toml > hybrid/inputs.sha256
"$JULIA" --project=env experiments/single_index/mave_hybrid_check.jl
"$JULIA" --heap-size-hint=1400M --project=env experiments/single_index/mave_hybrid_worker.jl 1
python3 - <<'PY_CHECK'
import csv
from pathlib import Path
h=list(csv.DictReader(Path('hybrid/part-0001.csv').open()));assert len(h)==16
assert all(r['status']=='ok' for r in h)
p=list(csv.DictReader(Path('confirmation/part-0001.csv').open()))
for rep in ['1','2']:
 a=next(r for r in h if r['replicate']==rep and r['method']=='stein_screen10')
 b=next(r for r in p if r['replicate']==rep and r['method']=='stein_screen10')
 assert abs(float(a['loss'])-float(b['loss']))<1e-10
 a=next(r for r in h if r['replicate']==rep and r['method']=='mave_select_rrr_refit')
 b=next(r for r in p if r['replicate']==rep and r['method']=='iterative_adaptive_lasso_refit')
 assert (a['missing'],a['extra'])==(b['missing'],b['extra'])
print('HYBRID_WORKER_CANARY_PASSED')
PY_CHECK
seq 2 4800 > hybrid/tasks.txt
parallel --jobs "$SLURM_NTASKS" --joblog hybrid/parallel.tsv --halt soon,fail=1 \
 srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project=env \
 experiments/single_index/mave_hybrid_worker.jl {} :::: hybrid/tasks.txt
python3 experiments/single_index/mave_hybrid_report.py
