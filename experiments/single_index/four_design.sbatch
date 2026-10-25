#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=4
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
STAGE=${1:?stage};PLAN="experiments/single_index/four_design_$STAGE.toml";OUT="$STAGE"
mkdir -p "$OUT"
sha256sum Project.toml Manifest.toml src/SingleIndexGraph/*.jl experiments/single_index/four_design* experiments/single_index/aggregation_common.jl experiments/single_index/multiprobe_common.jl > "$OUT/inputs-$SLURM_JOB_ID.sha256"
if [[ "$STAGE" == canary ]]; then
    "$JULIA" --project experiments/single_index/four_design_check.jl
fi
python3 - "$PLAN" "$OUT" > "$OUT/tasks-$SLURM_JOB_ID.txt" <<'PY'
import sys,tomllib,math
from pathlib import Path
plan=tomllib.loads(Path(sys.argv[1]).read_text())
for i in range(1,1+sum(math.ceil(s['replicates']/plan['shard_size']) for s in plan['scenario'])):
    if not (Path(sys.argv[2])/f'part-{i:04}.csv').exists():print(i)
PY
parallel --jobs "$SLURM_NTASKS" --joblog "$OUT/parallel-$SLURM_JOB_ID.tsv" --halt soon,fail=1 \
  srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project \
  experiments/single_index/four_design_worker.jl "$PLAN" {} "$OUT" :::: "$OUT/tasks-$SLURM_JOB_ID.txt"
python3 experiments/single_index/four_design_report.py "$PLAN" "$OUT"
