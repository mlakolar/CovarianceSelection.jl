#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/project2/mkolar_1314/SingleIndexGraph/20261024-sig-finalvalidation/.report-env/bin/python experiments/single_index/four_family_figures.py "$PWD"
sha256sum experiments/single_index/four_family_figures.py > reporting/figure_source.sha256
