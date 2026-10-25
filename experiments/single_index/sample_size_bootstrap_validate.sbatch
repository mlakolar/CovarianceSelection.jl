#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project -e 'include("experiments/single_index/sample_size_bootstrap.jl"); for task in 1:2; ssbootstrap("experiments/single_index/sample_size_bootstrap_canary.toml",task,"bootstrap_canary"); end'
python3 experiments/single_index/sample_size_bootstrap_assess.py experiments/single_index/sample_size_bootstrap_canary.toml bootstrap_canary
python3 test/sample_size_bootstrap.py bootstrap_canary
