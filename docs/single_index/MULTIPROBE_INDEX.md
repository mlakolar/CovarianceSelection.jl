# Multiple-probe Gaussian index study

This study estimates only the index. It does not fit a covariance curve or
alternate between a nonparametric fit and an index fit. Precision-matrix
coverage is outside its scope. The methods below are research extensions;
their empirical tuning rules are not established inference guarantees.

## Identification and estimators

Let Z have independent standard Gaussian coordinates and let
E[XX' | Z] = Sigma(Z'beta), with ||beta|| = 1. For a fixed response probe
Y_l = (q_l'X)(r_l'X), write m_l(t) = E[Y_l | Z'beta=t]. Then

```
M = E[Z Y'] = beta c',       c_l = E[m_l'(T)],
H_l = E[(ZZ' - I)Y_l] = b_l beta beta',
b_l = E[m_l''(T)].
```

The first-order matrix has rank one when c is nonzero. Individual probes may
be uninformative or have coefficients of opposite signs. Leading singular
directions pool squared signal rather than averaging signed columns. When
all first-order coefficients vanish, more first-order probes cannot repair
identification. Second-order identification separately requires b nonzero.

`probe_index_moments(Y, Z; order=1)` uses unbiased sample covariances between
probe responses and Stein features. Order two uses raw covariate products;
subtracting I disappears after sample response centering. The population
identity still requires the correct standard Gaussian covariate distribution.

`estimate_probe_index(Y, Z; method=..., support_size=...)` implements:

| Method | Construction |
|---|---|
| `:ols` | Multivariate OLS with intercept; leading left singular direction of slopes |
| `:moment` | Leading singular direction of first-order Stein moments |
| `:screen` | Rank full OLS slope rows by norm; retain k rows; refit all probe responses by OLS; take a singular direction |
| `:split` | Symmetric cross-product of independent first-order moment estimates |
| `:second` | Leading eigenvector of the sum of squared second-order moment matrices |
| `:second_split` | Leading eigenvector of symmetric independent-half second-order cross-products |

The default support size is d. No coefficient is required to be nonzero.
Directions use the largest-magnitude-coordinate sign convention, or a supplied
reference. Evaluation aligns signs with the truth only after estimation.

Probe responses are standardized to empirical unit variance. For direct
estimators all observations supply the scales. Standardized split estimators
reserve the first fifth of observations for calibration and split the remaining
rows in half. Conditional on independent calibration, cross-products have the
rank-one signal expectation. Reusing a moment half to learn a common scale
would destroy this exact conditional identity. Splitting therefore has a real
sample-size cost. With fixed external scales, use `standardize=false`; an
explicit two-way split must partition all supplied rows.

For second order the direct matrix is sum_l H_l H_l. For split estimates it is
sum_l (H_1l H_2l + H_2l H_1l)/2. Under independent halves and fixed scales,
its expectation is beta beta' sum_l b_l^2. Split estimates may be indefinite.
The estimator uses the largest positive eigenvalue, not the eigenvalue of
largest absolute magnitude. Invalid inputs, constant probes, and absent positive
signal eigenvalues produce explicit failures.

## Designs

Both designs use d=20 and Gaussian conditional responses. The primary beta has
four equal nonzero coefficients at coordinates 2:5. Thus its first coordinate
is truly zero. The mixed-sign sensitivity uses coefficients proportional to
(0.7, -0.5, 0.4, -0.3) at the same locations.

Write g(t) = tanh(t) + 0.5 tanh(t)^2. A signal multiplier a scales g.

* **Chain:** Omega_jj(t) = 1.5 + 0.15 a g(t), and adjacent edges are
  -(0.30 + 0.08 a g(t)). Other off-diagonal elements are zero.
  Diagonal dominance gives uniformly positive eigenvalues for all a tested.
* **Pairs:** Independent covariance blocks have unit diagonal and correlation
  rho(t) = 0.35 + 0.15 a g(t). All marginal variances are constant. Their inverse
  has fixed, nonzero pair edges. Square probes have no population signal.

There is no added privileged probe block. Probe construction never receives the
graph, the covariance curve, or beta. Pair locations are in the generator only.
The even-link diagnostic replaces g by tanh(t)^2. Gaussian quadrature checks
the first- and second-order identifying moment norms and positive definiteness.
These oracle calculations diagnose designs; they never tune an estimator.

## Probe dictionaries

* One randomly chosen coordinate product, as a baseline.
* All p coordinate squares.
* L random distinct coordinate products, L in {10,30,60}.
* L random projection squares with four random signed coefficients of magnitude
  1/2, for L in {10,30,60}.
* Mixed dictionaries: all p squares, L/2 products, L/2 projection squares.
* The complete quadratic dictionary: all p squares and all p(p-1)/2 products.

Probe seeds vary across repetitions and are independent of the data. Random
dictionaries are nested across L within a family. The complete dictionary is
an additional feasible benchmark for d=20. Its large response dimension can
still increase estimation noise and computation. It is evaluated with OLS,
first-order moments, and the two screening/refit choices, rather than expensive
second-order aggregation. Screening sizes are fixed at 6 and 10. All other
methods retain all 20 covariates. There are 75 fixed configurations.

## Stages and feasible tuning

1. A bounded Discovery canary checks moment identities, input failures,
   orientation, scale/shift invariance, representative runtime, and finite fits.
2. Development uses 300 paired repetitions for each combination of design and
   n in {100,200,400,800,1600,3200}, with p=40. Even-link diagnostics use 200
   repetitions per design at n in {400,800,1600}. Total: 4,800 datasets.
3. The initial development rule freezes six candidates: the two smallest average
   projective losses within each of the OLS, screening/refit, and spectral
   families. The twelve main scenarios receive equal weight. The spectral
   family contains moment, split, second, and second_split. Even-link cases
   cannot influence the shortlist. The overall smallest development average
   defines the fixed default.
   Development revealed a design-specific specialist that this global average
   omitted: complete quadratic OLS substantially outperformed random dictionaries
   for constant-variance pairs. Before inspecting any confirmation estimates,
   the initial confirmation was cancelled, its partial output retained, and the
   candidate pool expanded with complete OLS and its k=6 and k=10 refits.
   Final confirmation restarts every scenario with a fresh seed namespace.
   The nine-candidate pool and amended plan are frozen before this restart.
4. Confirmation uses fresh seeds and 1,000 repetitions per setting. It includes
   the twelve main settings, six p=160 settings (n=400,800,1600), eight signal
   sensitivities (a=0.5,1.5; n=400,1600; both designs), and four mixed-sign
   sensitivities (n=400,1600; both designs). Total: 30,000 datasets.

Confirmation compares every fixed shortlisted candidate, feasible validation
selection, and an explicitly labeled oracle choice among the same full-data
fits. The oracle is an evaluation benchmark and cannot be deployed.

Feasible tuning uses three folds. Fit candidates on two folds. Divide the
held-out fold into two subsets and score a candidate direction u by

```
u' M_validation1 M_validation2' u.
```

Every candidate is scored with the SAME complete quadratic reference dictionary,
whose scales are learned only on training observations. This prevents candidates
from improving their score merely by changing the validation dictionary or its
dimension. Conditional on training, the score expectation is
||c_reference||^2 (u'beta)^2. Scores may be negative; they are not clipped.
Average the three scores, choose the largest, and refit on all observations.
Failed validation candidates receive -Inf. If all fail, record a selector failure.
Truth never enters feasible selection. This rule is proposed and must earn its
place through independent confirmation, especially at n=100 or weak signal.

## Evaluation and execution

Primary loss: 1 - (beta_hat'beta)^2. Index prediction loss under a fresh standard
Gaussian covariate equals 2 - 2|beta_hat'beta|. Also retain angle, missing and
extra variables, runtime, and reasons for failures. Summaries include medians,
90th percentiles, and Monte Carlo standard errors. Failures remain in denominators
with maximal projective loss one. There is no exact-support gate.

Fixed-candidate runtimes measure fitting on cached probe responses; validation
runtime includes its own probe construction and fits. These two timings serve
different purposes and must not be compared as identical costs. Worker timings
can include initial JIT compilation and are descriptive, not controlled benchmarks.

All numerical fitting, integration, aggregation, and figures run on Discovery.
The launcher uses GNU Parallel within one Slurm allocation, with eight single-CPU
workers and 16 GB total RAM initially. The remaining final confirmation work
can resume with 32 single-CPU workers and 64 GB, under the user's existing larger
allocation authorization. Atomic finalized shards are reused without refitting;
only unfinished shards run. Scientific settings and seeds do not change.
BLAS and Julia threads are one per worker. A canary
uses one CPU and 4 GB. Plans, source hashes, task outcomes, per-repetition results,
and the lightweight run registry supply the recreation recipe. No routine data
archive or unrelated filesystem audit is required.

Entry points:

```
python3 experiments/single_index/multiprobe_plans.py development development.toml
sbatch experiments/single_index/multiprobe.sbatch canary development.toml canary
sbatch experiments/single_index/multiprobe.sbatch development development.toml development
python3 experiments/single_index/multiprobe_plans.py confirmation_specialists experiments/single_index/multiprobe_confirmation_final.toml development/shortlist.json
sbatch experiments/single_index/multiprobe.sbatch confirmation experiments/single_index/multiprobe_confirmation_final.toml confirmation_final
```

The source snapshot, installed Julia environment, and manifest used by a run
must be preserved for recreation. Results and job identities are added after
validation; scheduler completion alone is not scientific validation.
