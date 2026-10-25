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
module load ver/2506 gcc/14.3.0
python3 experiments/single_index/mave_report.py experiments/single_index/mave_development.toml development
python3 experiments/single_index/mave_report.py experiments/single_index/mave_confirmation.toml confirmation
python3 experiments/single_index/mave_interpret.py .
