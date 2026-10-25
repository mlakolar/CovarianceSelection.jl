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
for stage in known estimated; do
    if test -f "$stage/summary.csv"; then
        python3 experiments/single_index/sample_size_report.py "$stage"
    fi
done
/project2/mkolar_1314/SingleIndexGraph/20261024-sig-finalvalidation/.report-env/bin/python experiments/single_index/sample_size_plot.py .
python3 experiments/single_index/sample_size_followup_plans.py \
 experiments/single_index/sample_size_known.toml known/summary.csv experiments/single_index
