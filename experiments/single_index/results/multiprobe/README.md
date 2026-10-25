# Gaussian multiple-probe index experiments

Primary simulations are complete. All scientific computation and plotting ran
on Discovery. There are 4,800 development datasets (360,000 records) and 30,000
fresh confirmation datasets (330,000 records). Counts, identities, and numerical
failures were checked; both stages have zero failed fits.

The study estimates the index only. It uses d=20, conditionally Gaussian
responses, fixed precision support, a truly zero first index coefficient, and
no privileged informative probe block. There is no exact-support gate.

The full design and interpretation are in
[MULTIPROBE_INDEX.md](../../../../docs/single_index/MULTIPROBE_INDEX.md).

## Main comparison

Mean projective error at p=40,n=3200; smaller is better. Each design has 1,000
paired repetitions. An unrelated random direction has expected error 0.95.

| Procedure | Chain | Constant-variance pairs |
|---|---:|---:|
| Projection probes, OLS | 0.278 | 0.801 |
| Projection probes, select 10 and refit OLS | 0.217 | 0.794 |
| Complete quadratics, OLS | 0.832 | 0.189 |
| Complete quadratics, select 10 and refit OLS | 0.834 | 0.140 |
| Feasible held-out selection | 0.416 | 0.251 |
| Oracle choice among nine fitted candidates | 0.153 | 0.093 |

Probe weighting matters substantially. Liberal selection and OLS refitting
help with useful probes. Aggressive selection can improve the median while
worsening mean and upper-tail errors. The proposed held-out tuning rule is not
reliable enough to recommend. At p=160,n=1600, even oracle choice among the nine
candidates has mean error 0.764 for constant-variance pairs. Second-order Stein
aggregation did not solve the even-link diagnostic.

These outcomes do not establish a fundamental sample-size barrier. The complete
quadratic dictionary spans the projection probes; its poorer performance can
reflect inefficient aggregation of correlated noise. No method from this study
is promoted as a universal first stage for precision inference.

## Files

* `development/summary.csv`: all 75 configurations, including adverse outcomes.
* `development/shortlist.json` and `FROZEN*.txt`: development choices and freeze records.
* `confirmation_final/summary.csv`: all fixed candidates, feasible selection, and oracle choice.
* `confirmation_final/paired_controls.csv`: refit versus OLS on identical probes,
  and feasible selection versus every fixed candidate.
* `confirmation_final/selection_counts.csv`: feasible selected configurations.
* `confirmation_final/random_product_signal.csv`: exact probabilities that a
  random product dictionary misses every informative pair. At L=60 the
  probabilities are 0.198 for p=40 and 0.684 for p=160.
* `confirmation_final/resume_identity.json`: identical scientific fitting hashes
  across production attempts.
* `confirmation_final/recovery_identity.json`: 16 regenerated partial rows match
  their earlier values; only 54 unfinished datasets were regenerated.
* `figures/`: standalone PDF and PNG figures, with Monte Carlo error bars for
  mean losses. The gray dotted line marks random-direction error. The even-link
  figure uses development repetitions, not fresh confirmation.

## Selected-support diagnostic correction

Primary direction errors, prediction errors, angles, failures, and tuning choices
are validated. Secondary support-count columns are withheld from these versioned
tables pending selected-row replay. Dense eigensolvers left tiny values outside
the selected support, which made numerical-nonzero counts inappropriate.
The public estimator now zeros those entries exactly; the worker counts
membership in `fit.selected`.

Job `13424115` is queued for the diagnostic replay and new exact-zero tests.
Discovery currently reports `QOSMaxJobsPerUserLimit`. The replay checks normal
equations against direct QR and verifies that every primary summary quantity
remains exactly unchanged. It does not rerun primary direction estimation.
Original raw records, including their deprecated numerical-nonzero diagnostics,
remain intact on Discovery. Use the primary columns here until this correction
has completed and been reconciled into the reporting copy.

## Provenance

Primary root:
`/project2/mkolar_1314/SingleIndexGraph/20261025-sig-multiprobe`.
Diagnostic replay root:
`/project2/mkolar_1314/SingleIndexGraph/20261025-sig-multiprobe-support`.

Development job: `13423559`. Final confirmation: `13423599`, resource-only
resume `13423691`, and tail recovery `13423994`. Reporting and 138 focused
regression checks: `13423645`. Canary jobs: `13423525` and `13423553`.
The earlier writer canary `13423550` failed and was repaired before production.
The six-candidate initial confirmation `13423572` was cancelled before inspecting
its estimates; its 135 finalized partial shards are retained but excluded.
The nine-candidate confirmation uses an independent fresh seed namespace.

Source snapshots, manifests, raw result shards, and attempt logs remain under
the two roots above. The project run registry supplies launch commands and
resource history. This directory versions compact reporting inputs rather than
the full raw dataset records. Recreate stages with the supplied TOML plans and
Slurm launchers; do not pool the excluded initial confirmation with final results.
