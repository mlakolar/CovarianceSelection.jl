# Gaussian MAVE index comparison

All scientific computation ran on Discovery. Z has 20 standard Gaussian
coordinates. The response dimensions are 20 and 40. The four families and
n=200/400/800/1600 settings are described in
[the method and findings](../../../../docs/single_index/MAVE_REFINEMENT.md).
No precision estimation or inference is included.

* `development/`: 40 paired replications per setting and three bandwidth factors.
* `confirmation/`: 300 fresh paired replications per setting at factor one.
  The two failed holdout-arm records remain in `failures.csv` with loss one.
* `hybrid/`: the same confirmation datasets, with feasible support selection
  followed by Stein/RRR or OLS coefficient refitting. `oracle_support_rrr`
  receives true support and is a diagnostic only.
* `interpretation/`: complete tables, paired Monte Carlo errors, selection and
  trajectory consistency checks, and the visually inspected comparison figure.
* `convergence/`: 16 predetermined longer runs; all converged within 411 updates.
* `holdout_guard/`: replay of the single validation draw behind the two failures.
  One candidate had a nonfinite prediction. The guarded selector returns a
  finite two-variable fit with large error. These results do not replace the
  primary failures or repair its generally poor support choices.
* `smoke/` and `tiny/`: numerical/worker validation fixtures, not additional
  confirmation evidence. The smoke plan is the historical pre-extension plan.
* `environment/`: the exact experiment project and lockfile. Absolute local
  package paths are replaced when `Pkg.develop` stages a recreation environment.
* `run_manifest.json`: jobs, resource requests, source snapshots, and validation.

Primary loss is `1 - abs(beta' * estimate)^2`. Finite capped fits retain their
measured error. Failed fits receive loss one. There is no support gate.
Timing columns contain compilation and shared-stage effects and are not
controlled speed benchmarks.

Recreate on Discovery in an isolated run directory. Stage this repository's
sources and the committed KernelSmoothing snapshot there. Create its dependency
archive from a clone without using or modifying a dirty checkout:

```sh
git -C /path/to/KernelSmoothing.jl archive --format=tar \
  --output=/path/to/run/kernel-357e261.tar \
  357e26152026a3e1fb427acd221f5981c968722e Project.toml src
```

Submit `mave_setup.sbatch`, then `mave.sbatch development` and
`mave.sbatch confirmation`. `mave_hybrid.sbatch` validates its first shard
against the primary confirmation before launching its remaining workers.
`mave_interpret.sbatch` validates complete stage outputs and generates the
paired report. The plotting subprocess reuses the project's `.report-env`;
`SIG_PLOT_PYTHON` can select an equivalent Python with Matplotlib.
Use the Slurm profile and task-scoped SSH workflow from the project run record.
Numerical work and summaries must run under Slurm.

The historical primary source is `source-v3.tar.gz` in the recorded Discovery
run directory. Its pre-guard helper is also in commit `bd26779`. The hybrid
source snapshot and later report hashes are recorded separately. Current code
contains the finite-validation guard and public hybrid API; the retained primary
comparison predates the guard. Raw shards, full trajectories, interrupted
checkpoint files, and source snapshots remain on Discovery.

Committed CSV files use LF line endings. Remote CSV-writer artifacts may use
CRLF; field values are unchanged. Recorded transfer hashes refer to the
original remote bundles.
