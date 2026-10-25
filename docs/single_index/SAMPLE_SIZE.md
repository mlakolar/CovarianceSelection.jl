# Sample size and generated-index uncertainty

This follow-up separates precision inference with a known index from the cost
of estimating that index. It applies no index-support advancement gate. The
earlier requirement of 99% exact index-support recovery was a simulation-design
choice, not an established sample-size requirement for precision inference.

## Frozen initial design

Use the fixed-support Gaussian design in `GAUSSIAN_FIXED_SUPPORT.md`. Covariates
are Gaussian. Responses are Gaussian conditional on the covariates; their
marginal distribution need not be Gaussian. The true unit direction has four
components equal to 0.5. The informative probe remains prespecified. This study
does not evaluate selection of a random probe.

The known-index ladder crosses n = 100, 200, 400, 800, 1600, 3200 with p = 40,
160, 400. A second ladder uses p = 160 and disjoint hub blocks with at most six
leaves. All 24 settings use 1000 independent replications. Settings with the same
sample size share seeds across dimensions; their response arrays need not be
nested. Report locations -0.5, 0, 0.5 and a diagonal, active edge and zero edge.

For each hub block, set the diagonal to 1.5+0.1*tanh(t) and its active edges to
-(0.38+0.06*tanh(2t)). The support and signs are fixed for every real t. The
smallest main-block eigenvalue is at least 1.4-0.44*sqrt(6)>0. Hub targets are
(1,1), (1,2), (2,3); the last pair consists of two leaves. The probe is unchanged.
The hub sensitivity changes conditioning and diagonals as well as degree. It
does not identify a pure sparsity effect.

Use fixed normal-score penalty C=1.1, probability 0.1 per nuisance regression,
and the realized weighted predictor normalization. Use final bandwidth
2.5*n^(-1/5), pilot bandwidth 5*n^(-1/5), and local quadratic smoothing.
Compare union-support WLS refitting, the existing degrees-of-freedom/HC3/inverse
bias corrections, and true nuisance coefficients. All use the known index.
No penalty is selected by cross-validation.

## Reporting

Coverage denominators include every requested replication. Failed fits count
as noncoverage. Bias, empirical SD, SE and finite widths describe available
fits; report their failure counts alongside them. Report results separately by
location and entry. Wilson intervals quantify Monte Carlo uncertainty.

`effective_n` describes nonnegative pilot weights. `final_effective_n` equals
the inverse sum of squared local-polynomial intercept weights, which can be
negative individually. `final_local_n` counts observations inside the final
kernel window. These quantities have different roles.

Report the full coverage curves. Any observed sample-size transition is specific
to this graph, signal, bandwidth and target location. Investigate final and pilot
bandwidths separately near the transition, then confirm promising changes with
independent replications. Do not label a setting successful solely because its
intervals are wide or because support recovery is exact.

## Estimated-index comparison

At p=40, d=20 and n=200, 400, 800, 1600, 3200, compare the known direction,
the frozen ordinary-mean selector with multiplier 1.25 followed by OLS, and OLS
on the true index support. The last method is infeasible. Each setting has 1000
replications; methods share observations. Keep the precision estimator fixed
to the union-refitted rule. Do not exclude failures or condition advancement on
index support. This distinguishes direction-estimation error from variable
selection and preserves the small-sample failures.

## Candidate adjustment

The pairs bootstrap resamples complete (X,Z) observations. Its full version
repeats index support selection, OLS normalization, nuisance support selection,
union refitting and final smoothing. Its fixed-index version holds the original
estimated direction fixed but repeats the second stage on the same resamples.
Both use the original deterministic penalty and bandwidth rules.

For each entry and location form the studentized pivot
(estimate_star-estimate)/se_star. If its empirical 0.025 and 0.975 quantiles
are q_low and q_high, the bootstrap interval is
[estimate-q_high*se, estimate-q_low*se]. Compare with the original Wald interval.
Each mode uses the same resampled rows. With a known direction their outputs
must coincide, which supplies a direct software check.

Do not discard invalid bootstrap draws. A bootstrap interval is unavailable if
any requested draw has an invalid pivot for that entry/location. Count it as
noncoverage and report invalid draws and interval availability. This strict
policy exposes instability; it is not a remedy for failed first stages. Report
coverage among available intervals only as a secondary diagnostic.
The invalid-pivot count equals the requested draw count when the original fit
fails; resampling is then skipped. It does not count executed solver calls.

This is an exploratory candidate, not a proved confidence procedure. Index
error can cause bias as well as extra variability. Resampling after hard
selection may be unreliable near selection boundaries. Favorable coverage in
these Gaussian settings would not establish validity under weak identification
or general covariate distributions.

## Execution

All fitting, simulation summaries and scientific plots run under Slurm on
Discovery. Plans are generated by `sample_size_plans.py`. The launcher is
`sample_size.sbatch`; workers are `gaussian_fixed_worker.jl` and
`sample_size_bootstrap.jl`. Run records live in the parent project's lightweight
registry under `20261025-sig-sample-size`. A source snapshot and the reused Julia
1.13.2 manifest identify the actual inputs.

## Primary results

The known-index ladder completed 24,000 datasets and 648,000 entry/location/
method rows. The estimated-index comparison completed 5,000 datasets and
135,000 rows. Every expected shard and row was present. Identities were unique;
failure outcomes remained in coverage denominators.

For the main union-refitted estimator, ranges below are over three locations.
They are empirical proportions from 1000 replications, not confidence limits.

| Graph | n | p | Diagonal | Active edge | Zero edge |
|---|---:|---:|---:|---:|---:|
| Chain | 100 | 400 | 85.0--89.7% | 93.2--95.2% | 95.4--95.7% |
| Chain | 200 | 40 | 93.1--94.3% | 94.0--94.3% | 93.5--94.8% |
| Chain | 200 | 400 | 91.2--91.9% | 94.8--95.0% | 94.5--95.4% |
| Chain | 400 | 40 | 94.1--95.0% | 94.8--95.3% | 94.6--96.0% |
| Chain | 400 | 400 | 93.7--94.5% | 94.8--95.3% | 94.9--96.3% |
| Hub, degree 6 | 400 | 160 | 92.7--93.9% | 93.7--95.2% | 94.8--96.0% |
| Hub, degree 6 | 800 | 160 | 93.4--95.0% | 94.3--95.2% | 92.9--96.4% |
| Hub, degree 6 | 1600 | 160 | 95.0--96.1% | 94.3--95.3% | 94.3--94.9% |

The 400-response chain at n=400 has final effective size about 73--83 across
locations. The hub at n=1600 has about 225--255. The favorable chain needs far
fewer observations than the earlier index-support gate suggested. The hub
sensitivity needs more. There is no universal threshold, and observed coverage
is not monotone in n at every entry/location. Near 95%, the Monte Carlo standard
error for 1000 replications is about 0.7 percentage points. Mild deficits remain.

The frozen feasible index rule has failures in 978, 906, 651, 143 and 2 of 1000
replications at n=200,400,800,1600,3200 respectively. This failure means the
thresholded direction is zero or has a zero orientation anchor. These are not
failures of the precision solver. At n=1600 its unconditional pointwise coverage
is 80.5--82.0% across the nine targets. At n=3200 it is 93.7--95.3%.

OLS using the true index support has no such failures. At n=400 its diagonal
coverage is 93.4--94.9%, active-edge coverage 95.2--95.5%, and zero-edge coverage
94.0--95.3%. This infeasible comparison shows that variable selection and index
normalization can dominate the sample-size requirement in this design.
It does not establish a feasible procedure at that sample size.
The frozen selector was chosen for exact support recovery in larger samples.
Its poor small-sample behavior is not a lower bound on the sample size needed
to estimate an index. Less aggressive screening and refitting with a few extra
variables remain distinct possibilities; this frozen comparison does not test
them.

Primary figures and full reporting inputs are in
`experiments/single_index/results/sample_size/{known,estimated}/`.

## Independent follow-up specifications

For each graph/dimension family, locate the first n where the observed minimum
coverage of the primary estimator is at least 93% and remains at least 93% at
every larger tested n. This is a development locator, not a declaration of
95% validity. Investigate that n and its predecessor. It selects n=100,200 for
the p=40 chain; n=200,400 for the p=160 and p=400 chains; n=800,1600 for the hub.
The exact minima and chosen sizes are saved in
`sample_size_followup_selection.json`.

At each chosen size, use 500 independent replications and final/pilot factor
pairs (0.7,2), (1,2), (1.3,2), (1,1), (1,3), relative to 2.5*n^(-1/5).
Thus each change modifies one bandwidth only. Compare all three original
variants. Seeds begin at 947000000; scenario IDs separate specifications.
These comparisons are sensitivity analysis, not an independently confirmed
coverage-optimal tuning rule.

The sensitivity run completed 20,000 independent datasets and 540,000 result
rows in 400 shards. In the p=400 chain at n=400, the fresh baseline gives
93.8--95.4% diagonal coverage. Reducing the pilot factor from 2 to 1 gives
89.8--91.8%; the minimum nuisance-support containment rate falls from 77.0%
to 21.0%. Increasing the pilot factor to 3 gives 92.6--95.2% coverage and
95.2% minimum containment. Better containment alone does not guarantee better
coverage. The scenario seeds differ across specifications; these are independent
comparisons, with Monte Carlo standard errors about one percentage point.

At n=800 in the hub design, pilot factor 1 gives 89.4--94.4% diagonal coverage
and 11.2% minimum containment. Factors 2 and 3 give 93.6--94.2% and 93.6--95.2%
coverage, with minimum containment 90.4% and 99.6%. The pilot needs enough local
information for nuisance selection and refitting. These results support retaining
the wider pilot as a benchmark; they do not establish its optimality. Narrowing
the final smoother alone does not consistently repair diagonal undercoverage.

Very small windows also produce occasional nonpositive residual-covariance
failures, including oracle-nuisance fits. Those outcomes remain in the tables.
The correction bundle is not uniformly preferable to the primary estimator.

The bootstrap uses 300 independent outer replications at n=800,1600,3200,
399 pairs draws, and the feasible frozen index and true-support OLS diagnostic.
Seeds begin at 967000000. Outer Monte Carlo uncertainty is about 1.3 percentage
points near 95%. The experiment is exploratory; modest coverage differences
require additional independent confirmation.

## Bootstrap findings

All 900 outer datasets completed, producing 48,600 rows in 90 shards. The
following ranges are across nine entry/location pairs. Availability is the
proportion of outer datasets where an interval can be constructed under the
strict invalid-draw policy above.

| n | Feasible-index Wald coverage | Full-bootstrap availability | Full-bootstrap coverage |
|---:|---:|---:|---:|
| 800 | 34.7--38.0% | 0% | 0% |
| 1600 | 81.3--84.7% | 2.0% | 1.3--2.0% |
| 3200 | 92.0--95.7% | 43.0% | 38.0--41.3% |

At n=3200 the original feasible index is available in 99.3% of these datasets.
Nevertheless, repeating hard selection in 399 draws often produces at least
one undefined direction. Rejecting an entire interval in that case explains
much of the full bootstrap's poor unconditional coverage. These results do not
show that every bootstrap treatment of selection must fail; they describe this
frozen estimator and explicitly stated failure policy.

Instability from undefined indices is not the only issue. With OLS on the true
index support, all intervals are available. Full-bootstrap coverage is
86.7--92.7%, 89.7--93.7%, and 89.7--94.7% at n=800,1600,3200. Corresponding Wald
ranges are 92.7--97.3%, 93.0--96.7%, and 92.3--96.3%. The full bootstrap does not
improve this diagnostic. Merely re-estimating both stages in a pairs bootstrap
is therefore not a validated adjustment here. Fixed-index bootstrap-t intervals
also show no consistent improvement over Wald intervals.

The next methodological issue is to separate first-stage bias, selection
instability and the studentized bootstrap's calibration. This study does not
resolve it. Paired coverage changes and their Monte Carlo standard errors are
saved in `bootstrap/paired_coverage.csv`; do not choose an adjustment from the
most favorable individual cell. No bootstrap change is promoted to a default.

The numerical validation passed 92 existing inference checks, 107 focused
design/penalty/worker checks, 26 hub checks, and paired-bootstrap identity checks
at 18 dataset/entry/location combinations. A plotting job failed on a tiny
negative error-bar length caused by roundoff at zero coverage. Only plotting
was repaired; simulation outputs and coverage summaries were preserved.
