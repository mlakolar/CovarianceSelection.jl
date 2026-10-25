# Fixed-support Gaussian study

Status: index comparison, pointwise development/confirmation, multiplier bands,
and paired grid refinement completed and validated on Discovery.
This study supersedes the earlier sign-changing block generator for this
experiment only. Existing studies and software defaults are unchanged.

## Design

The covariates are independent standard Gaussian. The unit index direction has
four components equal to 0.5 and all other components zero. One sensitivity uses
components 1/sqrt(2), 1/sqrt(6), 1/sqrt(6), 1/sqrt(6). Responses conditional on
covariates are mean-zero Gaussian, with covariance the inverse of the precision
specified below. Marginal responses need not be Gaussian.

Partition the first p-2 responses into chains of four, with a final shorter
chain when necessary. Set their precision diagonals to 1.25+0.10*tanh(t) and
chain edges to -(0.38+0.06*tanh(2*t)). All other entries are zero. Active edge
magnitudes lie in [0.32,0.44], diagonals in [1.15,1.35], and the smallest
main-block eigenvalue is at least 0.27. The graph and edge signs are constant
for every real index. The generator uses each chain adjacency eigensystem to
sample its covariance without a full p-dimensional factorization per row.

The last two responses form an independent conditional probe block with unit
variances and correlation rho(t)=0.55+0.40*tanh(2*t). Its precision is
[1 -rho; -rho 1]/(1-rho^2). The edge never vanishes. Its response product times
the known Gaussian covariate score identifies the index direction. Probe
identification is separate from the graph entries used for inference.

## Index comparison and advancement gate

Compare ordinary-mean and five-block median-of-means moments, each with soft
thresholding, hard thresholding, or support selection followed by OLS of the
probe product on the selected covariates and an intercept. All directions are
unit-normalized and oriented by the first coordinate. OLS on the true support
is an explicitly infeasible diagnostic. Threshold multipliers are
[0.75,1,1.25,1.5,2] times the existing empirical coordinate-scale rate.

Use 1000 replications in each setting:

| n | p | d | Direction |
|---:|---:|---:|---|
| 1200 | 40 | 20 | equal |
| 2400 | 40 | 20 | equal |
| 2400 | 80 | 50 | equal |
| 2400 | 40 | 20 | unequal |
| 4800 | 40 | 20 | equal |
| 4800 | 80 | 50 | equal |
| 9600 | 80 | 50 | equal |

The larger settings are a prespecified sample-size ladder for locating a
favorable regime. They are not substituted silently for lower-sample failures.
Record squared direction error, exact support, containment, missed and extra
variables, runtime and failures. A setting advances only with observed exact
support >=99%, containment >=99.5%, and no numerical failures for the selected
feasible candidate. Wilson intervals accompany these empirical proportions;
the gate is not a proof of a population probability. Prefer the feasible candidate passing the largest number of equal-direction
settings, breaking ties by mean squared direction error across its passing
settings. Advance at most two settings: its smallest and largest passing
sample-size regimes. Retain the unequal-direction case as a sensitivity,
without allowing it to choose the primary rule. Freeze the selected method and threshold
before independent second-stage replications. Never select an index using its
true direction separately within an inference replication.

Index development seeds: 811000000 + scenario_id*100000 + replication.
Canary seeds are in the disjoint 809000000 namespace.

## Fixed penalties and subsequent stages

Use refitted common-complement square-root lasso through the existing scaled-
lasso solver. Normalized pilot weights are w_i. Predictor loadings remain
sqrt(sum_i w_i X_ij^2). For m predictors and tail probability alpha=0.1, set

    lambda = C * Phi_inverse(1-alpha/(2*m))
               * max_j sqrt(sum_i w_i^2 X_ij^2 / sum_i w_i X_ij^2).

This is a fixed safety constant and a realized design scale, with no penalty CV.
It incorporates the Kish effective size and predictor-specific squared-weight
variation. It approximates homoskedastic normalized-score calibration; it does
not remove local model approximation or estimated-index error. The optional
penalty_tests argument divides the tail budget across a prespecified collection
of fits. Its default is one per regression, not a uniform nuisance guarantee.

After the index gate passes, compare C=[0.75,1,1.1,1.25], pilot factors
[1.5,2,3], and final factors [0.6,1,1.4], multiplying 2.5*n^(-1/5).
No lambda CV is used. Assess three points [-0.5,0,0.5] and entries (1,1), (1,2),
(1,4). The zero entry is within a chain, so it generally has nonzero marginal
covariance. Compare separate-column support refits, common-union refits, and
union refits with residual degrees-of-freedom, HC3 and inverse-curvature
corrections. All endpoints are additive Wald endpoints. Known indices and true
nuisance coefficients are diagnostic controls.

Pointwise development uses independent seeds and 200 replications. Freeze a
small number of penalty/bandwidth specifications before 1000-replication
confirmation. Keep coverage separate by point and entry and report bias,
empirical SD, RMS estimated SE, widths, failures and nuisance support diagnostics.
Do not advance an index setting that fails the gate.

Multiplier bands follow pointwise assessment, using 1999 shared Gaussian
observation multipliers and a 21-point grid on [-0.75,0.75], with a smaller
paired 41-point refinement. Curvewise and joint coverage are distinct.
Finite-grid performance does not establish continuum coverage.

## Completed index comparison

Discovery job 13422307 completed all 350 shards, 7000 datasets and 217000
candidate rows. Counts, unique identities and failure retention passed.
Job 13422310 applied the index-only advancement rule. Software validation
passed 92 existing inference checks and 103 new design/penalty/worker checks.

The frozen feasible candidate is ordinary-mean support selection with threshold
multiplier 1.25, followed by OLS on selected covariates and an intercept. The
threshold is a development choice; it is not reselected using truth within
coverage replications. Lower-sample failures remain in the comparison.

| Setting | Exact support | Contains all true variables | Mean squared direction error |
|---|---:|---:|---:|
| n=1200, p=40, d=20 | 14.8% | 14.9% | 0.34187* |
| n=2400, p=40, d=20 | 91.2% | 92.2% | 0.03733* |
| n=2400, p=80, d=50 | 81.2% | 81.2% | 0.06061* |
| n=4800, p=40, d=20 | 99.1% | 100% | 0.01056 |
| n=4800, p=80, d=50 | 99.6% | 100% | 0.01025 |
| n=9600, p=80, d=50 | 99.1% | 100% | 0.00553 |

Asterisks denote error averages over successful fits only. Support proportions
include failures. The first three cases had 353, 12 and 42 zero-direction or
anchor failures respectively. They were not advanced. At n=2400 with unequal
components, exact recovery was 59.9%; this sensitivity also was not advanced.

At n=4800,p=40,d=20, direction MSE was 0.04008 for soft thresholding,
0.01335 for hard thresholding, 0.01056 after OLS refitting, and 0.01021 for
OLS with the true support. This comparison uses the same observations and the
same threshold. The post-selection fit is close to its oracle-support benchmark
in this favorable setting. This does not establish a general support theorem.

The 99.1% exact-support proportion has Wilson interval approximately
[98.30%,99.53%]; observed 99% is not proof that the population probability
exceeds 99%. All missed-support, extra-variable and failure outcomes are saved.
The selected rule advances n=4800,p=40,d=20 and n=9600,p=80,d=50, the smallest
and largest passing sample-size regimes. Pointwise development uses independent
829000000 seeds and job 13422320.

Compact index reporting inputs are in
`experiments/single_index/results/gaussian_fixed_index/`. Raw shards and source
snapshots remain in the task-owned Discovery run directories recorded in the
run registry.

## Central confirmation specification

Before seeing the independent 853000000 outcomes, freeze the central pilot
factor 2 and final factor 1. Compare fixed C values 0.75, 1.1 and 1.25 and all
three refit/correction variants, with 1000 replications per passing setting.
This central specification is fixed while development is running; no empirical
coverage optimum is selected from its 36-point tuning surface. The development
grid supplies bandwidth and penalty sensitivity, and confirmation assesses the
central rule. Keep actual index-support performance in the new replications.
Uniform-band production follows pointwise assessment, with its outcomes below.

## Independent pointwise confirmation

Jobs 13422352 and 13422356 completed 1000 replications per setting, including
paired true-index/true-nuisance controls. The main study retained 324000 rows
in 100 shards; the oracle control retained 18000 rows in 40 shards. Job
13422357 computed summaries. Job 13422367 passed 107 focused regression,
design, penalty and worker checks after the diagnostics repair.

At C=1.1, pilot factor 2, final factor 1, ordinary-mean/OLS index estimation,
and common-union WLS refitting without extra corrections, coverage ranges over
three evaluation points were:

| Setting | Diagonal | Active edge | Zero edge |
|---|---:|---:|---:|
| n=4800,p=40,d=20 | 94.2–95.1% | 93.3–94.7% | 94.3–95.7% |
| n=9600,p=80,d=50 | 93.3–94.9% | 93.9–94.3% | 94.6–95.2% |

These are pointwise ranges, not simultaneous coverage. With the true index,
the corresponding diagonal ranges were 94.9–95.2% and 94.0–95.8%. The
true-nuisance/true-index diagonal ranges were 94.8–95.3% and 93.9–95.7%.
Thus mild deviations persist even for oracle nuisances; neither selection nor
lambda alone explains every deviation. At 1000 replications, a coverage near
95% has Monte Carlo standard error about 0.7 percentage points. Wilson
intervals and each individual point remain in summary.csv. No scenario or
entry is dropped because its coverage is unfavorable.

The estimated-index exact-support proportions were 99.1% and 99.4% in this
independent confirmation. For central C=1.1, nuisance support containment was
100% across the recorded point/entry fits, with mean extra union sizes 0.063
and 0.050. Adding residual-maker, HC3 and inverse-curvature corrections changed
coverage little: estimated-index diagonal ranges became 94.1–94.9% and
93.5–94.9%. They are not required to obtain the principal favorable-design
behavior, and are not shown to solve every remaining deficit.

The 200-replication development surface compares all four constants and both
bandwidth grids. At the central bandwidths, C=0.75 selected about 1.45 and 1.80
extra union variables per fit; C=1.1 selected about 0.048 in both settings;
C=1.25 selected about 0.008 and 0.005. Coverage and widths changed much less
than these support sizes. This supports fixed score calibration followed by
refitting as a useful starting point in this design. It is not a universal
constant or bandwidth theorem. Confirmation comparisons of C=0.75,1.1,1.25
are retained separately.

Reporting inputs are in `results/gaussian_fixed_development/` and
`results/gaussian_fixed_confirmation/` beneath `experiments/single_index/`.
The 877000000 band namespace is independent of index development, pointwise
development and pointwise confirmation. Its primary configuration preserves
C=1.1, pilot factor 2 and final factor 1. Union-refitted and corrected procedures
are compared to known-index and genuine true-nuisance controls on a 21-point
grid with 1999 multipliers. A paired 41-point check uses the first 200 datasets
of the smaller sample setting. Bands do not retroactively certify nominal
pointwise coverage or continuum coverage.

The pilot-factor range 1.5–3 had little effect on the central estimated-index
coverage surface or mean widths. The final factor changed widths substantially:
at n=4800, mean width was 0.3084, 0.2388 and 0.2023 for factors 0.6, 1 and 1.4
at C=1.1 and pilot factor 2. Smaller final bandwidths did not uniformly improve
coverage. In the 200-replication n=9600 development case, the smallest final
factor had one diagonal coverage of 90.0%; it is retained in the surface.
These experiments do not support an automatic aggressive-undersmoothing rule.

## Independent multiplier bands

Jobs 13422382 and 13422384 completed the band study and paired grid refinement;
13422387 reported saved results. The main study has 700000 rows in 400 shards,
including genuine known-index/true-nuisance controls. The 41-point check has
26000 rows in 20 shards. Counts and duplicate-identity checks passed. Every
requested band fit succeeded; no unsuccessful fit was removed.

For ordinary-mean/OLS index estimation, common-union refitting and no extra
finite-sample correction, nominal 95% coverages on 21 grid points were:

| Setting | Diagonal curve | Active-edge curve | Zero-edge curve | Joint family |
|---|---:|---:|---:|---:|
| n=4800,p=40,d=20 | 95.0% | 95.4% | 95.0% | 94.7% |
| n=9600,p=80,d=50 | 94.1% | 95.3% | 96.4% | 95.4% |

Curvewise columns use separate critical values. The family column uses one
larger critical value across every entry and grid point. It is not the
intersection of the three curvewise bands. Joint Wilson intervals are
[93.13%,95.93%] and [93.92%,96.53%]. Corresponding true-index family coverage
was 95.2% and 95.1%; true-nuisance/true-index family coverage was also 95.2%
and 95.1%. With the extra finite-sample corrections, feasible family coverage
was 94.6% and 95.5%. These small changes do not establish a correction benefit
in this favorable design.

Exact index support in the independent band namespace was 98.8% and 99.3%.
All true variables were retained in every replication; the discrepancies were
extra variables, with means 0.012 and 0.007. This illustrates why an observed
99% development gate must not be interpreted as a guaranteed population
probability. No threshold was changed after these outcomes.

On the first 200 paired smaller-setting datasets, joint coverage changed from
92.5% on 21 points to 93.0% on 41 points: difference +0.5 percentage points,
paired Monte Carlo SE 0.5 points. Diagonal and active-edge curve coverages were
unchanged; zero-edge curve coverage changed from 94.0% to 93.0% (paired SE
0.71 points). Mean family width changed from 0.3927 to 0.3946. The entire
paired comparison is retained. This is a resolution diagnostic, not a
continuum certificate.

The practical starting specification supported by these experiments is an
ordinary-mean, support-selected OLS index refit; fixed normal-score C=1.1;
common-union nuisance WLS refits; pilot bandwidth 5*n^(-1/5); and final
bandwidth 2.5*n^(-1/5), on the standardized Gaussian index scale. It remains
a design-specific empirical choice. Pointwise coverage retains mild deficits,
and no general tuning, weak-signal or continuum theorem is asserted.

All numerical work ran on Discovery. Two initial technical attempts failed
(missing dependency depot; variable scope in new support diagnostics); both
were repaired before production. Their logs and source versions are retained.
All final simulation and reporting jobs completed with exit code 0:0. The
module defaults and earlier simulation results remain unchanged.

Local CSV reporting copies use LF line endings. Their scientific cell values
are unchanged from the Discovery exports.
