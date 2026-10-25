#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=4
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load ver/2506 gcc/14.3.0 parallel/20250522
# All workers and numerical summaries remain inside this one Slurm allocation.
# Each stage has immutable finalized shards and resumes only missing shards.
for stage in grid paired corrections confirmation double_training; do
    if [[ "$stage" == corrections && ! -f experiments/single_index/promising_coverage_corrections.toml ]]; then
        continue
    fi
    printf '%s %s START\n' "$(date -u +%FT%TZ)" "$stage" >> stages.tsv
    bash experiments/single_index/promising_coverage.sbatch "$stage"
    printf '%s %s VALIDATED\n' "$(date -u +%FT%TZ)" "$stage" >> stages.tsv
    if [[ "$stage" == grid ]]; then
        python3 experiments/single_index/promising_coverage_select.py grid > selection.log
    fi
done
printf '%s COMPLETE\n' "$(date -u +%FT%TZ)" >> stages.tsv
