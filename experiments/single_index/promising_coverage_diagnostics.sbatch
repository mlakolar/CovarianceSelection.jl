#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 experiments/single_index/promising_coverage_interpret.py "$PWD"
python3 experiments/single_index/promising_coverage_orientation.py "$PWD"
sha256sum experiments/single_index/promising_coverage_interpret.py experiments/single_index/promising_coverage_orientation.py > reporting/diagnostic-source.sha256
