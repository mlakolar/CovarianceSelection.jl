# Iterative index refinement and variable selection

This is an index-estimation experiment. It does not fit precision matrices or
construct confidence intervals. The implementation is a rank-one, projected
multivariate MAVE adaptation for quadratic response probes. Its empirical
performance and numerical convergence are separate questions. No new statistical
or support-recovery theorem is asserted.

The Discovery comparison is complete. Adaptive MAVE is useful for covariate
selection, followed by a first-order Stein/RRR coefficient refit. One dense
update changes little, while dense iteration usually increases direction error.
A simpler fixed-threshold Stein refit remains competitive in two designs.
See the [comparison figure](../../experiments/single_index/results/mave/interpretation/index-comparison.png)
and [complete paired tables](../../experiments/single_index/results/mave/interpretation/comparisons.md).

## Model, observations, and probe banks

Write `Z` for the n-by-20 Gaussian covariate matrix and `X` for the n-by-p
conditionally Gaussian response matrix. Rows are independent. The direction
`beta` has unit norm, with entries 2--5 equal to `(1,1,-1,-1)/2` and all other
entries zero. In particular, its first coordinate is zero. The index is
`T = Z * beta`. Conditional response means are zero, and conditional covariance
is `Sigma(T)`. Response coordinates are independently permuted in each dataset.
No estimator sees the permutation, true support, covariance graph, derivative,
or true direction.

Use the four families already specified in [FOUR_GAUSSIAN_DESIGNS.md](FOUR_GAUSSIAN_DESIGNS.md).
Here are their strongest-amplitude versions. Let `A` be the path adjacency
matrix divided by two. Let `B` be a fixed perfect matching without path edges,
chosen to make `A*B != B*A`. Both have operator norm at most one.

| Family | Conditional covariance or precision |
| --- | --- |
| Correlation | Independent 2-by-2 covariance blocks with diagonal one and correlation `0.15 + 0.55*tanh(T)` |
| Opposed scale | Independent blocks `exp(s*0.6*tanh(T))*[1 .25; .25 1]`, with balanced alternating block signs `s` |
| Fixed precision diagonal | `Omega(T) = I + (0.1 + 0.6*tanh(T))*A` |
| Noncommuting | `Omega(T) = I + 0.2*A + 0.6*tanh(T)*B` |

The latter two families have constant precision diagonals. Edges can cross zero.
These are the existing informative designs, rather than a new fixed-support
restriction. In every case the conditional covariance is positive definite.

Use p=20 and 40, and n=200, 400, 800, and 1600. Only the index-training sample
is generated; n is the total sample size in this study. Use the same observations
and probe dictionary for the paired methods. The opposed-scale dictionary is
all coordinate squares, 30 randomly chosen distinct coordinate products, and
30 random squared projections onto four response coordinates. The other designs
use every coordinate square and distinct coordinate product. These bank choices
were fixed from the earlier independent Stein study. They are design-specific;
this experiment does not establish a universally selectable bank.

Let `Y` collect these scalar probes. Standardize each column by its sample mean
and sample standard deviation once, obtaining `V`. All iterations retain these
scales. There is no response-covariance inverse weighting. With Gaussian rows of
Z, the first-order Stein identity aligns the cross moment between Z and Y with
the index, when the vector of average link derivatives is nonzero.

## Shared local regressions and direction update

Initialize with the existing screened/refitted Stein reduced-rank regression:
rank rows of the multivariate OLS slope matrix, keep the ten largest row norms,
refit OLS on those variables, and extract the whitened leading rank-one direction.
A dense 20-variable Stein fit is a second baseline. Screening ten variables is
an initializer, not an eligibility gate.

For current unit direction b, compute `t = Z*b`. The bandwidth is

$$h = c_h n^{-1/5}\operatorname{sd}(t).$$

The main factor is `c_h=1`. Development also evaluates 0.7 and 1.5. At each
iteration discard the bottom and top 5% of projected observations as regression
anchors. All observations, including those in the tails, remain neighbors.
For anchor i, normalize Gaussian kernel weights over j and solve

$$\min_{a_i,s_i}\sum_j w_{ij}\|V_j-a_i-s_i(t_j-t_i)\|_2^2.$$

The intercept and slope are vectors of probe means and derivatives. One pair of
smoothing matrices supplies all probe fits. `KernelSmoothing.smoothermatrix`
(the name at the frozen dependency commit) is an independent implementation of
these matrices. It is also used for validation-sample prediction. The core MAVE
routine fuses the same degree-one Gaussian formulas with the neighborhood
calculation. Numerical equivalence is checked on all four families.

Freeze neighborhoods, intercepts, and slopes. Solve the quadratic problem

$$\min_u \frac{1}{2mL}\sum_{i,j}w_{ij}
  \|V_j-a_i-s_i u^T(Z_j-Z_i)\|_2^2,$$

where m is the number of anchors and L the number of probes. The normal matrix
and right-hand side are

$$Q=\frac{1}{mL}\sum_{i,j}w_{ij}\|s_i\|_2^2
       (Z_j-Z_i)(Z_j-Z_i)^T,$$

$$r=\frac{1}{mL}\sum_{i,j}w_{ij}(Z_j-Z_i)s_i^T(V_j-a_i).$$

A relative diagonal ridge of `1e-10 * trace(Q)/20` stabilizes the solve.
Set the next direction to `u / norm(u)`, aligning its sign to the previous
iterate. This is an unconstrained update followed by normalization. Multiplying
the frozen local slopes by `norm(u)` preserves its fitted products. It is not a
claim that normalization solves a quadratic problem with fixed slopes on the
unit sphere. Refit the local curves before the next update. Because neighborhoods
change, the full alternating objective need not decrease between iterations.
The frozen quadratic objective before and after every raw update is recorded.

Compare exactly one update with iteration to a sign-aligned direction change
below `1e-4`, with a cap of 200 updates for dense MAVE and 60 for the penalized and restricted
fits. Dense trajectories retain every iterate to distinguish slow convergence
from deterioration. An iterate at the cap is labeled
`maxiter`, retained, and reported separately from convergence. All 20 variables
are updated in the dense procedure, including variables omitted at initialization.
The sign convention does not identify the model's intrinsic sign ambiguity.
Direction errors are therefore sign invariant.

## Variable selection

Start from the dense MAVE fit. Explore three fixed selection rules, followed by
unpenalized restricted MAVE, with the same bandwidth and stopping rule.

* **Threshold/refit:** retain `abs(b[j]) > C*sqrt(log(20)/n)`, with
  C=0.5, 1, and 1.5 in development and C=1 in confirmation.
* **Lasso/refit:** minimize `u'Q*u/2-r'u + lambda*sum(abs(u))` at the dense fit.
  Set `lambda = f * max(abs(r))`.
* **Adaptive Lasso/refit:** weight the absolute penalty by
  `v[j] = 1/max(abs(b[j]),0.05)`, and use
  `lambda = f * max(abs(r[j])/v[j])`. The positive floor permits any coordinate
  to enter. For both Lasso rules, f=0.03, 0.1, and 0.25 in development;
  confirmation uses the prespecified f=0.1.

These levels are exploratory, not Belloni--Chernozhukov calibrated penalties.
Their coordinate scales assume comparable covariate units; here Z is standard
Gaussian. The bounded adaptive weights do not establish an oracle property.
They do not use cross-validation or the true support. Coordinate descent checks
the penalized quadratic's KKT residual. An empty model is a recorded failure;
there is no oracle fallback. Restricted refitting removes shrinkage but can
still retain a wrong model. It cannot recover a variable excluded by the final
selection, so omissions must be measured explicitly.

In addition, apply ordinary or adaptive Lasso at **every** direction update,
starting from screened Stein. Adaptive weights are fixed from that initializer
with the same 0.05 floor. Recompute `lambda_max` at each update, use fraction
0.1, and keep every covariate eligible. After stopping, refit the final support
without a penalty. Record penalized-stage and final-refit convergence separately.
This separates early regularization from selection after a potentially
overfit dense iteration.

Also explore an honest validation split. Use 75% of observations to initialize,
fit dense MAVE, rank all twenty coefficients, and refit nested supports of
sizes 2, 4, 6, 10, and 20. Compare standardized probe prediction error on the
other 25%, with all centering, scaling, feature rankings, and curves fitted on
the training part. Choose the minimum validation error, or separately the
smallest support within one standard error of that minimum. Refit the selected
support on all n observations. This is a feasible support-size selector, not
an oracle comparison. The usual one-standard-error rule may be too aggressive
because irrelevant probes dilute validation differences. Retain this possibility
as an adverse outcome.

## Refitting index coefficients after selecting variables

The default sparse procedure refits the nonparametric MAVE criterion without a
penalty. Development suggested separating support selection from coefficient
estimation. Accordingly, a follow-up uses MAVE to select covariates, then
estimates their index coefficients using the first-order Stein relation.

For fixed probes, write `g(t) = E[Y | T=t]`. With standard Gaussian covariates,

$$E[ZY^T]=\beta E[g'(T)]^T.$$

Thus a multivariate OLS slope matrix has population rank one when the vector
of average derivatives is nonzero. The same identity holds after fixed probe
rescaling. Sample centering and scaling are plug-in steps; no exact finite-sample
unbiasedness claim is made for sample-standardized responses.

On a data-selected covariate support, use the following calculations. `Zs` is
the centered selected design, and `V` is the sample-standardized probe matrix.

```julia
G = Zs' * Zs / (n-1)
M = Zs' * V / (n-1)
L = cholesky(Symmetric(G)).L
N = L \ M
# OLS: the leading eigenvector of (L' \ N) * (L' \ N)'.
# RRR: the leading eigenvector u of N*N', then normalize(L' \ u).
```

Embed the resulting direction into the full 20 coordinates. Both refits are
unpenalized. OLS and RRR use all fixed probes; neither chooses an informative
probe from the generating covariance graph. Align to the initializer when
possible. An orthogonal new support uses the largest-component sign convention
instead of producing an orientation error.

The follow-up has 300 replications in all 32 confirmation settings, using the
same observations as the primary comparison. The feasible methods are:

* adaptive MAVE selection for 60 updates at fraction 0.1, followed by OLS or RRR;
* one adaptive penalized MAVE update at fraction 0.1, followed by RRR;
* thresholding the Stein initializer at C=0.5, 1, or 1.5 times
  `sqrt(log(20)/n)`, followed by RRR.

All choices were frozen from development before their confirmation summaries.
Selection and refitting reuse the n index-training observations. A separately
labeled true-support RRR fit is an infeasible diagnostic. It never supplies
support to a feasible method or an eligibility gate. The simpler thresholding
arm tests whether iterative MAVE adds value beyond inexpensive Stein selection.
The final paired report checks that the two coefficient-refit arms have the
same MAVE selection outcomes on matched replications.

## Comparisons, failures, and reproducibility

Primary loss is squared sine direction error `1 - (beta'bhat)^2`. Report mean
and Monte Carlo standard error, median, 90th and 95th percentiles, the fraction
above 0.25, and paired differences from screened Stein. Also report the
sign-aligned index prediction loss, exact-support frequency, all-true-variable
containment, average omissions and extras, update count, final change, and time.
A numerical failure receives loss one in aggregate direction comparisons.
Finite `maxiter` fits retain their measured loss and are counted as nonconverged.
There is no support gate and no selective removal of difficult replicates.
Timing fields include some per-worker compilation and omit data/probe generation.
They are computational diagnostics, not controlled solver benchmarks.

The bounded canary checks the normal matrix and right-hand side against a
separately accumulated stacked regression, Lasso KKT conditions, noiseless
recovery, sign equivariance, restricted support, reentry of an initially missing
variable, a truly zero first coordinate, and local smoother equivalence with
KernelSmoothing. Simulation smoke tests precede development and confirmation.
Development and confirmation use disjoint seed streams. Confirmation levels
are prespecified; no per-replicate truth-based tuning is allowed.

Run under Discovery Slurm using `mave_setup.sbatch`, followed by
`mave.sbatch smoke`, `mave.sbatch development`, and `mave.sbatch confirmation`.
The experiment environment has local path dependencies on this frozen source
and KernelSmoothing commit `357e26152026a3e1fb427acd221f5981c968722e`.
The actively edited KernelSmoothing checkout is not used or modified.
The lightweight project registry is
`SingleIndexGraph/.slurm/runs/20261026-sig-mave-refinement.json`.
Workers preserve finalized shards and resume only missing shards.

The relevant multivariate reference is
[Zhang (2021), weighted MAVE and adaptive group Lasso](https://doi.org/10.1016/j.jmva.2021.104753).
For a single index, row-group selection reduces to a weighted absolute penalty.
Our errors are correlated and heteroscedastic quadratic probes. Published
asymptotic and selection results are motivation, not established guarantees
for this implementation or its tuning rules.

## Execution and results

All scientific computation ran under Discovery Slurm, with Julia 1.13.2.
The updated canary (job 13447860) passed 105 core, smoothing, and regression
checks and all 44 expected worker fits. Public hybrid options passed 64 checks,
including agreement with the frozen experiment. The hybrid worker passed nine
helper checks and matched the primary first-two-replication canary.

Development job 13447940 validated 1280 datasets and 49920 fits, with zero
numerical failures. It retained 15729 finite capped fits. Confirmation uses
9600 fresh datasets: 300 replications in each of 32 settings. Its 105600 fits
include two failed holdout-selector arms, both retained with loss one, and
30233 finite capped fits. The hybrid study reuses those same 9600 datasets;
it has 76800 reported fits and no numerical failures. It is not another 9600
independent datasets.

Confirmation job 13447950 was checkpointed after 1576 finalized shards. Only
missing shards resumed in disjoint allocations 13448632 and 13448636, which
completed in 35m57s and 34m47s. The hybrid job 13448768 completed in 26m12s.
Its Slurm transport warnings did not prevent any worker or parent-job
completion. Expected identities, row counts, and unique replicate/method keys
passed. All failed-arm records remain available. No support gate was used.

Reporting jobs 13448654 and 13449137 passed 19200 paired consistency checks.
The frozen MAVE selection outcomes agree across the MAVE-coefficient and
Stein-coefficient refit arms. Every stored dense trajectory matches its final
fit: 491831 development and 1292151 confirmation trajectory records. The
comparison figure was visually inspected after rendering.

At p=40, d=20, n=800, the confirmation mean squared-sine losses are:

| Design | Stein10 | One dense update | Dense MAVE, cap 200 | Adaptive selection + MAVE refit | Adaptive selection + RRR refit | Stein threshold C=1.5 + RRR |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Correlation | 0.01517 | 0.01537 | 0.09834 | 0.01809 | 0.00479 | 0.00680 |
| Fixed precision diagonal | 0.01735 | 0.01741 | 0.07830 | 0.01488 | 0.00532 | 0.00743 |
| Noncommuting | 0.00807 | 0.00816 | 0.03604 | 0.00542 | 0.00249 | 0.00261 |
| Opposed scale | 0.00810 | 0.00803 | 0.01825 | 0.00426 | 0.00344 | 0.00254 |

Thus coefficient refitting matters even after the same support is selected.
RRR is preferable to the nonparametric MAVE coefficient refit in these cases.
It also improves on OLS coefficient refitting in correlation, fixed diagonal,
and noncommuting; OLS and RRR are essentially tied in opposed scale. The
first-order derivative condition is essential to the Stein/RRR hybrid.

At these same settings, the adaptive MAVE selector retained every true variable
in all 300 replications for each design:

| Design | Exact support | Average extra variables | Penalized stage at its 60-update cap |
| --- | ---: | ---: | ---: |
| Correlation | 98.7% | 0.013 | 280/300 |
| Fixed precision diagonal | 99.0% | 0.010 | 234/300 |
| Noncommuting | 100.0% | 0.000 | 55/300 |
| Opposed scale | 78.7% | 0.240 | 0/300 |

These are empirical strong-design results for a fixed 60-update selector.
They are not convergence or support-recovery guarantees. Finite penalized
iterates often have a useful support before their coefficient changes meet
the numerical tolerance. The final RRR refit is a closed-form calculation.

The fixed C=1.5 threshold baseline is important. At n=800, its RRR refit is
better than the MAVE hybrid in opposed scale and close in noncommuting.
The MAVE hybrid improves on it in correlation and fixed diagonal by paired
loss differences -0.00200 (MCSE 0.000275) and -0.00211 (MCSE 0.000264).
All C values, including 0.5 and 1, remain in the tables. These constants are
experimental choices, not universally calibrated penalties.

One dense update produces little improvement over Stein. In development,
increasing the bandwidth factor from one to 1.5 reduces dense MAVE error and
improves convergence, but it still exceeds the Stein10 mean loss in all four
p=40,n=800 cases. A predetermined longer check used the first two development
replications per design at p=40,n=400/800. All 16 fits converged within 411
updates; 14 ended with larger direction error than their Stein initializers.
This separates statistical deterioration from merely hitting the primary cap.
It does not establish a general failure of other MAVE implementations or weights.

The particular holdout probe-prediction selector is poor here. At p=40,n=800,
minimum-validation-loss mean direction errors range from 0.077 to 0.194.
The usual one-standard-error rule often drops two true variables and has mean
losses about 0.43--0.46. Prediction-risk tuning of the whole noisy probe bank
should not be presumed to select an accurate index.

The two numerical failures came from one correlation draw at p=40,n=1600,
replication 179. A single nonfinite prediction invalidated the four-variable
validation model. The old selector propagated its NaN score into an invalid
index operation. A finite-candidate guard now excludes unavailable candidates
and gives an explicit error if none are available. Three guard checks passed.
Replaying the exact draw chose a two-variable model with loss 0.5041 for both
selection rules. The guard repairs availability, not its poor variable choice.
The primary tables retain the original two failed-arm records; the repaired
case is a separately labeled diagnostic in `results/mave/holdout_guard/`.
The historical pre-guard helper is in commit `bd26779` and the frozen run
snapshot. Other primary index procedures are unchanged by this guard.

Small-sample limitations remain. At p=40,n=200, the RRR hybrid's mean loss is
0.514 in correlation and 0.349 in fixed diagonal. At n=400 it is 0.055 and
0.039. Reliable support at n=800 does not establish that n=200 is adequate,
or that this procedure works for weak indices, other covariance paths, or
high-dimensional covariates. Only d=20 is studied here. No precision estimation
or inference was run in this follow-up.

Compact results, paired Monte Carlo errors, failures, upper tails, provenance,
and the figure are retained under
[`experiments/single_index/results/mave/`](../../experiments/single_index/results/mave/).
Raw shards and trajectories remain on Discovery. The plotting environment is
the existing project `.report-env`, with one numerical-library thread; its
version is recorded alongside the figure. The run registry records submission,
checkpoint, and source-snapshot details. Publication sources also provide the
new public hybrid options; the frozen primary simulation source predates them.

## Public API

```julia
using CovarianceSelection.SingleIndexGraph

# Y is a fixed matrix of quadratic probes; Z has one covariate vector per row.
# This public example uses the OLS screening initializer. The study uses RRR.
b0 = estimate_probe_index(Y, Z; method=:screen, support_size=10).direction
one = refine_probe_index(Y, Z, b0; updates=1)
full = refine_probe_index(Y, Z, b0; updates=200)
sparse = select_mave_variables(Y, Z, b0; method=:adaptive_lasso,
                               level=0.1, iterative=true, refit_updates=60,
                               refit_method=:rrr)
# sparse.selection_status describes penalized iteration.
# sparse.status describes the final unpenalized support refit.
```

A new independent dataset would be projected as `Znew * sparse.direction`.
Use `refit_method=:rrr` or `:ols` in `select_mave_variables` for the coefficient
refits described above. The default `:mave` retains the nonparametric refit.
The output direction is normalized but is identified only up to sign. This API
does not choose an orientation for precision inference.
