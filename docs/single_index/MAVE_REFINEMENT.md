# Iterative index refinement and variable selection

This is an index-estimation experiment. It does not fit precision matrices or
construct confidence intervals. The implementation is a rank-one, projected
multivariate MAVE adaptation for quadratic response probes. Its empirical
performance and numerical convergence are separate questions. No new statistical
or support-recovery theorem is asserted.

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

## Comparisons, failures, and reproducibility

Primary loss is squared sine direction error `1 - (beta'bhat)^2`. Report mean
and Monte Carlo standard error, median, 90th and 95th percentiles, the fraction
above 0.25, and paired differences from screened Stein. Also report the
sign-aligned index prediction loss, exact-support frequency, all-true-variable
containment, average omissions and extras, update count, final change, and time.
A numerical failure receives loss one in aggregate direction comparisons.
Finite `maxiter` fits retain their measured loss and are counted as nonconverged.
There is no support gate and no selective removal of difficult replicates.

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

The core and KernelSmoothing canary passed 51 checks on Discovery, job 13447709.
The initial 32-dataset smoke run (job 13447770) produced all 288 expected
fits, with no numerical failures. It retained 95 finite fits at the 60-update
cap. The dense cap was therefore extended to 200 for the main comparison.
The updated canary (job 13447860) passed 105 numerical and regression checks
and all 44 expected worker records. Its four datasets are validation, not
substantive simulation evidence.

The development plan has 40 paired replications in each of 32 settings and
three bandwidth factors. The independent confirmation plan has 300 paired
replications per setting and uses the prespecified factor one. Both compare
Stein, one update, dense MAVE, and the feasible variable-selection procedures.
The fixed confirmation levels were not selected from its true index errors.
Jobs 13447940 and 13447950 are running on Discovery, with 128 single-core
workers each. Full-node requests were cancelled while pending and replaced
with smaller per-node footprints; no numerical outputs were discarded.
No inference is authorized in this campaign.
