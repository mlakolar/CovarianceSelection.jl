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
module load ver/2506 gcc/14.3.0
export JULIA_DEPOT_PATH="$PWD/depot:/project2/mkolar_1314/SingleIndexGraph/20261025-sig-index-confirmation/depot:/home1/mkolar/juliaPkg/depot:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
mkdir -p env dependencies/KernelSmoothing
# Committed dependency snapshot, rather than the owner's migration in progress.
tar -xf kernel-357e261.tar -C dependencies/KernelSmoothing
if [[ ! -f env/Project.toml ]]; then
 if [[ -f experiments/single_index/results/mave/environment/Manifest.toml ]]; then
  cp experiments/single_index/results/mave/environment/Manifest.toml env/Manifest.toml
 else
  cp /project2/mkolar_1314/SingleIndexGraph/20261026-sig-promising-coverage/Manifest.toml env/Manifest.toml
 fi
 cat > env/Project.toml <<'TOML'
[deps]
CovarianceSelection = "fb9f034c-9bf8-4e47-a0a4-723273363029"
TOML
fi
"$JULIA" --project=env -e 'using Pkg; Pkg.develop([PackageSpec(path=pwd()),PackageSpec(path=joinpath(pwd(),"dependencies","KernelSmoothing"))]); Pkg.instantiate()'
"$JULIA" --project=env experiments/single_index/mave_check.jl
