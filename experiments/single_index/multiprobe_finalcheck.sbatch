#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
sha256sum experiments/single_index/multiprobe_figures.py test/probe_index.jl \
 experiments/single_index/multiprobe_finalcheck.sbatch src/SingleIndexGraph/*.jl \
 > "report-source-$SLURM_JOB_ID.sha256"
"$JULIA" --project -e 'using Test,CovarianceSelection; include("test/single_index.jl"); include("test/index_orientation.jl"); include("test/probe_index.jl")'
/project2/mkolar_1314/SingleIndexGraph/20261024-sig-finalvalidation/.report-env/bin/python \
  experiments/single_index/multiprobe_figures.py .
python3 - <<'PY'
from pathlib import Path
import shutil,json
out=Path('reporting');out.mkdir(exist_ok=True)
for folder in ('development','confirmation_final','canary','specialist_canary','smoke2'):
    source=Path(folder);destination=out/folder;destination.mkdir(exist_ok=True)
    for path in source.iterdir():
        if path.is_file() and not path.name.startswith('part-'):
            shutil.copy2(path,destination/path.name)
shutil.copytree('figures',out/'figures',dirs_exist_ok=True)
for job in ('13423525','13423550','13423553','13423559','13423572','13423599'):
    path=Path(f'slurm-{job}.out')
    if path.exists():shutil.copy2(path,out/path.name)
for path in Path('.').glob('*sha256'):shutil.copy2(path,out/path.name)
receipt=dict(initial_confirmation_partial_shards=len(list(Path('confirmation').glob('part-*.csv'))),
             initial_confirmation_location='confirmation',
             initial_confirmation_used_in_final_evidence=False,
             reason='Development-only specialist amendment before inspecting confirmation estimates.',
             raw_production_locations=['development','confirmation_final'],
             raw_production_retained=True)
(out/'retention.json').write_text(json.dumps(receipt,indent=2)+'\n')
PY
