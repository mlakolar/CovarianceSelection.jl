#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
MODE=${1:?mode}; PLAN=${2:?plan}; OUT=${3:?output}
mkdir -p "$OUT"
sha256sum Project.toml Manifest.toml src/SingleIndexGraph/*.jl \
 experiments/single_index/multiprobe* test/probe_index.jl "$PLAN" > "$OUT/inputs-$SLURM_JOB_ID.sha256"
if [[ "$MODE" == canary ]]; then
    "$JULIA" --project experiments/single_index/multiprobe_canary.jl "$OUT"
else
    python3 - "$PLAN" > "$OUT/tasks.txt" <<'PY'
import sys,math,tomllib
with open(sys.argv[1],'rb') as f: plan=tomllib.load(f)
for i in range(1,1+sum(math.ceil(c['replicates']/plan['shard_size']) for c in plan['scenario'])):
    print(i)
PY
    parallel --jobs "$SLURM_NTASKS" --joblog "$OUT/parallel-$SLURM_JOB_ID.tsv" --halt soon,fail=1 \
      srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project \
      experiments/single_index/multiprobe_worker.jl "$PLAN" {} "$OUT" :::: "$OUT/tasks.txt"
    python3 experiments/single_index/multiprobe_report.py "$PLAN" "$OUT"
fi
