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
if test "${1:?stage}" = bootstrap; then
    python3 experiments/single_index/sample_size_bootstrap_report.py bootstrap
else
    python3 experiments/single_index/sample_size_report.py bandwidth
fi
/project2/mkolar_1314/SingleIndexGraph/20261024-sig-finalvalidation/.report-env/bin/python experiments/single_index/sample_size_plot.py .
