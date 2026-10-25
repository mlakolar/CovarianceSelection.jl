#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PARENT=${1:?primary run}
module load ver/2506 gcc/14.3.0 parallel/20250522
export JULIA_DEPOT_PATH="$PWD/depot:$PARENT/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
sha256sum Project.toml Manifest.toml src/SingleIndexGraph/*.jl experiments/single_index/multiprobe* > inputs.sha256
"$JULIA" --project -e 'using Test,CovarianceSelection; include("test/single_index.jl"); include("test/index_orientation.jl"); include("test/probe_index.jl")'
python3 - "$PARENT" > tasks.tsv <<'PY'
import sys,tomllib,math
from pathlib import Path
parent=Path(sys.argv[1])
for stage in ('development','confirmation_final'):
    plan=tomllib.loads((parent/f'experiments/single_index/multiprobe_{stage}.toml').read_text())
    for i in range(1,1+sum(math.ceil(s['replicates']/plan['shard_size']) for s in plan['scenario'])):
        print(stage,i,sep='\t')
PY
parallel --colsep '\t' --jobs "$SLURM_NTASKS" --joblog "parallel-$SLURM_JOB_ID.tsv" --halt soon,fail=1 \
  srun --exclusive --exact -N1 -n1 -c1 --mem=1900M "$JULIA" --heap-size-hint=1400M --project \
  experiments/single_index/multiprobe_support_worker.jl \
  "$PARENT/experiments/single_index/multiprobe_{1}.toml" {2} "results/{1}" :::: tasks.tsv
for stage in development confirmation_final; do
    PLAN="$PARENT/experiments/single_index/multiprobe_$stage.toml"
    python3 experiments/single_index/multiprobe_support_assess.py "$PLAN" "results/$stage"
    cp "$PARENT/$stage/summary.csv" "before-$stage.csv"
    python3 experiments/single_index/multiprobe_report.py "$PLAN" "$PARENT/$stage" "$PWD/results/$stage"
done
python3 experiments/single_index/multiprobe_support_finalize.py . "$PARENT"
/project2/mkolar_1314/SingleIndexGraph/20261024-sig-finalvalidation/.report-env/bin/python \
  experiments/single_index/multiprobe_figures.py "$PARENT"
python3 - "$PARENT" <<'PY'
import sys,shutil
from pathlib import Path
parent=Path(sys.argv[1]);dest=parent/'reporting'
for stage in ('development','confirmation_final'):
    for path in (parent/stage).iterdir():
        if path.is_file() and not path.name.startswith('part-'):shutil.copy2(path,dest/stage/path.name)
    shutil.copy2(Path('results')/stage/'receipt.json',dest/stage/'support_receipt.json')
shutil.copytree(parent/'figures',dest/'figures',dirs_exist_ok=True)
shutil.copy2('correction_identity.json',dest/'support_correction_identity.json')
shutil.copy2('inputs.sha256',dest/'support_inputs.sha256')
PY
