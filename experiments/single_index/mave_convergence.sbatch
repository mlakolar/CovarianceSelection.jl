#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
mkdir -p convergence
seq 1 16 > convergence/tasks.txt
parallel --jobs 8 --joblog convergence/parallel.tsv --halt soon,fail=1 \
 srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project=env \
 experiments/single_index/mave_convergence.jl {} :::: convergence/tasks.txt
python3 - <<'PY'
from pathlib import Path
import csv,json
paths=sorted(Path('convergence').glob('part-*.csv'));assert len(paths)==16
rows=[r for p in paths for r in csv.DictReader(p.open())];assert len(rows)==16
with Path('convergence/summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(json.dumps(dict(cases=len(rows),converged=sum(r['status']=='converged' for r in rows),maxiter=sum(r['status']=='maxiter' for r in rows))))
PY
