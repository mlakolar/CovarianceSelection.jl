#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
sha256sum Project.toml Manifest.toml src/SingleIndexGraph/*.jl experiments/single_index/four_design* experiments/single_index/aggregation_common.jl experiments/single_index/multiprobe_common.jl > population-inputs.sha256
/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project experiments/single_index/four_design_population.jl population
