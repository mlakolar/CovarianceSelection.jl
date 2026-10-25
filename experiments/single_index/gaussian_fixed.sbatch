#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=2
#SBATCH --ntasks=32
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
PLAN=${1:?plan}; OUT=${2:?output}; mkdir -p "$OUT"
sha256sum Project.toml Manifest.toml src/SingleIndexGraph/*.jl experiments/single_index/gaussian_fixed_worker.jl experiments/single_index/dgp.jl "$PLAN" > "$OUT/inputs-$SLURM_JOB_ID.sha256"
"$JULIA" --project -e 'using TOML; p=TOML.parsefile(ARGS[1]); println(sum(cld(c["replicates"],p["shard_size"]) for c in p["scenario"]))' "$PLAN" > "$OUT/task-count.txt"
seq 1 "$(cat "$OUT/task-count.txt")" > "$OUT/tasks.txt"
parallel --jobs "$SLURM_NTASKS" --joblog "$OUT/parallel-$SLURM_JOB_ID.tsv" --halt soon,fail=1 \
 srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --project experiments/single_index/gaussian_fixed_worker.jl "$PLAN" {} "$OUT" :::: "$OUT/tasks.txt"
python3 experiments/single_index/gaussian_fixed_assess.py "$PLAN" "$OUT"
