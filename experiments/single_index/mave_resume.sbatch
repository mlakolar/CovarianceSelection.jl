#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=8
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# Disjoint parity shards permit two fragmented allocations to share a resumed
# confirmation stage. Workers and scientific inputs are unchanged.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
PARITY=${1:?parity zero or one};[[ "$PARITY" == 0 || "$PARITY" == 1 ]]
PLAN=experiments/single_index/mave_confirmation.toml;OUT=confirmation
sha256sum Project.toml env/Project.toml env/Manifest.toml src/SingleIndexGraph/*.jl experiments/single_index/mave_worker.jl experiments/single_index/mave_common.jl "$PLAN" > "$OUT/resume-inputs-$SLURM_JOB_ID.sha256"
python3 - "$PLAN" "$OUT" "$PARITY" > "$OUT/tasks-$SLURM_JOB_ID.txt" <<'PY'
import sys,tomllib,math
from pathlib import Path
p=tomllib.loads(Path(sys.argv[1]).read_text());parity=int(sys.argv[3])
for i in range(1,1+sum(math.ceil(s['replicates']/p['shard_size']) for s in p['scenario'])):
 if i%2==parity and not (Path(sys.argv[2])/f'part-{i:04}.csv').exists():print(i)
PY
parallel --jobs "$SLURM_NTASKS" --joblog "$OUT/parallel-$SLURM_JOB_ID.tsv" --halt soon,fail=1 \
 srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project=env \
 experiments/single_index/mave_worker.jl "$PLAN" {} "$OUT" :::: "$OUT/tasks-$SLURM_JOB_ID.txt"
# A separate dependency job validates the union after both parities finish.
