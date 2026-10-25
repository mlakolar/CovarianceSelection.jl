#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
PLAN=experiments/single_index/multiprobe_confirmation_final.toml
OUT=confirmation_final
RECOVERY=tail_recovery
mkdir -p "$RECOVERY/generated"
python3 experiments/single_index/multiprobe_recover_check.py
sha256sum Project.toml Manifest.toml src/SingleIndexGraph/*.jl \
 experiments/single_index/multiprobe_common.jl experiments/single_index/multiprobe_worker.jl \
 experiments/single_index/multiprobe_recover* "$PLAN" > "$OUT/inputs-$SLURM_JOB_ID.sha256"
python3 experiments/single_index/multiprobe_recover.py prepare "$PLAN" "$OUT" "$RECOVERY" "$SLURM_NTASKS"
if test -s "$RECOVERY/workers.txt"; then
    parallel --jobs "$SLURM_NTASKS" --joblog "$RECOVERY/parallel-$SLURM_JOB_ID.tsv" --halt soon,fail=1 \
      srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project \
      experiments/single_index/multiprobe_recover_worker.jl "$RECOVERY/one_rep.toml" \
      "$RECOVERY/worker-{}.txt" "$RECOVERY/generated" :::: "$RECOVERY/workers.txt"
fi
python3 experiments/single_index/multiprobe_recover.py merge "$PLAN" "$OUT" "$RECOVERY"
python3 experiments/single_index/multiprobe_report.py "$PLAN" "$OUT"
