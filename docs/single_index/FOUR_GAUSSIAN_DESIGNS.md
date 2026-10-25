# Four Gaussian covariance designs: index estimation only

## Scope and frozen comparison

The user approved these four designs and Discovery execution after population
checks on 2026-10-25. Precision estimation and inference are excluded. Existing
non-iterative aggregation estimators are reused without changing their public API.
No exact-support gate, truth-based tuning, or first-coordinate anchor is used.

The grid has four families, three amplitudes, response dimensions p=20,40,
index dimensions d=4,20, and sample sizes n=100,200,400,800,1600: 240 settings.
Each setting has 300 paired repetitions, for 72,000 datasets. All fixed estimator
rules are evaluated; there is no development-selected shortlist. New random seeds
are disjoint from earlier campaigns. An eight-dataset canary uses its own seeds.
A later independent confirmation would be appropriate before paper-level claims
or selecting a precision-inference procedure.

Z is standard Gaussian. The unit index has four coefficients (1,1,-1,-1)/2,
on coordinates 1:4 when d=4 and 2:5 when d=20. The first coefficient is truly
zero when d=20. X conditional on Z is centered Gaussian with covariance Sigma(t),
t=Z'beta. A fresh independent response-coordinate permutation hides the graph
in each dataset. Every feasible method receives the same data and probe bank;
the permutation and population link are not passed to estimators.

## Covariance families

### Correlation changes, constant marginal variances

Independent two-by-two covariance blocks have diagonal one and correlation
rho(t)=0.15+a tanh(t), a=0.15,0.35,0.55. Across the entire amplitude grid, covariance
eigenvalues lie in [0.30,1.70]. At the weakest level edges are nonzero at finite t;
at higher levels covariance and precision edges cross zero. Squares have no
population first-order signal. Products carry first-order signal. Second-order
Stein moments vanish by symmetry; this is a diagnostic limitation, not a failure
of implementation.

### Opposing variance changes

Independent covariance blocks are exp(s_b a tanh(t)) R, where R has diagonal one
and off-diagonal 0.25, alternating block signs s_b=+1,-1, and a=0.2,0.4,0.6.
Both p values have equally many positive and negative blocks. Covariance eigenvalues
are between 0.75 exp(-0.6) and 1.25 exp(0.6). Precision support remains constant.
Diagonal and edge precision functions both vary. The trace probe has zero average
first-order signal because opposing blocks cancel, although individual probes
identify the index. This tests aggregation of responses with opposite loading signs.

### Varying precision edges, fixed precision diagonal

Let A be the path adjacency matrix divided by two, so its operator norm is below
one. Define Omega(t)=I+(0.10+a tanh(t))A, a=0.2,0.4,0.6, and Sigma=Omega^{-1}.
Precision eigenvalues lie in [0.30,1.70]. The diagonal stays one; allowed edges
can change sign. This removes the simultaneous precision-diagonal variation
that weakened covariance-edge derivatives in the earlier chain design.

### Noncommuting sparse precision matrices

Use the same A and a disjoint perfect-matching adjacency B with operator norm
one. The matching is generated with a fixed graph seed, excluding path edges,
and is checked to have AB != BA. Define
Omega(t)=I+0.20A+a tanh(t)B, a=0.2,0.4,0.6. Precision eigenvalues lie in [0.20,1.80].
The union graph has maximum degree three. Matching edges cross zero at t=0.
Covariance matrices do not share a fixed eigenbasis. Response permutations prevent
estimators from being given the matching labels.

## Population checks before fitting

For all 24 family/amplitude/p combinations, compute:

- Global covariance eigenvalue extrema and condition bounds, including t=+/-infinity.
- Average first derivative E[Sigma'(T)], its Frobenius norm, and its norm after
  scaling by the inverse square root of E[Sigma(T)].
- The second Stein matrix E[(T^2-1)Sigma(T)].
- Mean Gaussian information 0.5 E tr[Omega(T)Sigma'(T)Omega(T)Sigma'(T)].
- Trace-probe first signal and the first signal of the known score probe.
- Agreement of E[Sigma'(T)] with E[T Sigma(T)] and finite-difference derivatives.

Population integration uses 1,201 normal-weighted nodes on [-7,7]. All families
must have positive eigenvalues, nonzero first-order matrix signal, and positive
information. Trace cancellation is checked explicitly in the first two families.
These checks do not rank procedures using the true index and do not constitute
finite-sample accuracy or coverage guarantees.

## Fixed estimators and controls

All three dictionaries are evaluated: 60 sparse projection squares; a mixed bank
of all coordinate squares, 30 random products, and 30 projection squares; and the
complete quadratic dictionary. No bank is given graph information.

For each bank, compare moment aggregation, multivariate OLS/SVD, reduced-rank
regression (RRR), covariance-weighted OLS/RRR at rho=0.05,0.2,0.5, a two-fold
cross-fitted learned combined probe at rho=0.2, and independent split-moment
aggregation. RRR accounts for empirical covariate covariance before extracting
its leading direction. The mixed bank also evaluates second-order Stein
aggregation. Direct regression methods retain all d coordinates or screen to 10
when d=20 and refit; d=4 uses all four. Weighting rules share the same selected
row set, avoiding a change of support as an explanation for weighting gains.

Two feasible scalar baselines use sum_j X_j^2 and sum_{j<k} X_j X_k, with OLS
and the same row-screening rule. The first is an intentional negative control
in families with zero trace signal. Their weak population signals are retained.

Infeasible controls use q(X)=-X'Omega'(0)X with scalar-response OLS, on either all
d coordinates or the true four-variable support. They are diagnostic controls,
never inputs to feasible fitting, tuning, or support selection. At d=4 these two
controls coincide because all index variables are active.

## Outcomes and execution

Report mean squared-sine direction loss, Monte Carlo standard error, median,
90th percentile, sign-aligned prediction error, selected-support diagnostics,
failure counts, and paired differences against matching unweighted baselines.
Failures receive direction loss one and remain in all comparisons.

On 1,000 independent Gaussian evaluation covariates, compare normalized kernel
weights at index zero using h=n^{-1/5} and multiples 0.75,1.5. Record displacement
relative to bandwidth. These measure neighborhood distortion, not precision bias
or interval coverage. Shared preprocessing/weight-eigendecomposition costs are
recorded separately from per-method fits; timings are not full implementation
benchmarks.

Run identity: 20261025-sig-four-designs. Frozen plans are four_design_canary.toml
and four_design_comparison.toml. Julia 1.13.2 and the existing Manifest/depot are
reused. The population job precedes the canary. Only validated population checks
and a successful canary permit production. Production requests 64 single-threaded
CPUs on four nodes, 128 GB total RAM, and a two-hour limit. GNU Parallel launches
exclusive srun workers. Atomic 25-dataset shards allow resuming missing work
without dropping numerical failures. Expected raw output is approximately 1--3 GB.
Keep raw shards on Discovery and version compact reporting inputs and receipts.
