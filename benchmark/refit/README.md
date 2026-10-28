# Precision refit benchmark

This benchmark compares the four `covsel_refit` methods for Gaussian maximum
likelihood on a known graph. It measures how graph structure, dimension, and
sample size affect the choice of method. CliqueTrees.jl supplies an outside
reference: maximum-determinant covariance completion by Newton's method on a
chordal embedding. Its precision estimate is projected onto the original graph.

## Cases

The full suite uses p = 100, 500, 1000, 2000 and n = p/2, 2p, with five graphs:

- `band`: edges at distances 1 through 3; this graph is chordal.
- `grid`: 10×10, 20×25, 25×40, and 40×50 grids.
- `er`: Erdős–Rényi graphs with expected degree 4.
- `glasso-sparse` and `glasso-dense`: nonzero supports of QUIC fits, using
  `solve!(X, LogDetLoss(S), NormL1(λ, offdiag), ProxNewton())` from the identity,
  with `offdiag = ones(p, p) - I`. The penalties are respectively
  `sqrt(log(p)/n)` and `0.5sqrt(log(p)/n)`. Supports are used even when QUIC
  reports nonconvergence. The edge-list cache includes QUIC's penalty,
  iterations, and convergence flag.

For each case, the data-generating precision matrix has an independent
Erdős–Rényi support with expected degree 4. Off-diagonal entries have random
signs and magnitudes uniform on (0.2, 0.4). Each diagonal equals its absolute
off-diagonal column sum plus 0.2. Diagonal dominance bounds its smallest
eigenvalue below by 0.2 before rescaling. The population covariance is then
scaled to unit diagonal; this rescaling preserves the precision support.
Generate n independent Gaussian observations and use `S = X'X/n`, without
centering or sample standardization.

The case seed is `20261028 + 100000p + 10n + graph_index`, where the graph
indices are the order above. Fixed ER fitting graphs use a separate stream
with seed `case_seed(p, p, "er") + 1`, so they are shared across sample sizes.
Every case can be recreated independently. Data generation and graph fitting
are outside all refit timings.

## Methods and accuracy

`ChordalCompletion()` runs only when `is_chordal(E)`. `DualCoordinate()` and
`PrimalCoordinate()` use their defaults (`tol = 1e-10`, `maxiter = 1000`).
Exceptions, including dual breakdowns with singular sample covariances, are
recorded without retrying that method on the case. A returned nonconverged fit
keeps status `ok` and its separate convergence flag.

For p ≤ 500, an untimed primal fit sets the ADMM target to
`max(primal_kkt, 1e-9)`. Untimed ADMM fits try equal absolute and relative
tolerances 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, in that order, with 20,000 iterations.
Choose the first tolerance meeting the target. If none does, use 1e-10 and
mark the calibration unmatched. A failed or nonfinite primal certificate also
leaves calibration unmatched; the fallback target is 1e-9. Timed fits start
fresh with the chosen settings.

At p ≥ 1000, reuse the tolerance from p = 500 for the same graph and n/p
ratio, and time one ADMM fit. The iteration cap is 600 at p = 1000 and 300
at p = 2000. Each iteration costs an eigendecomposition. The p = 2000
budget is therefore at most 300 eigendecompositions, potentially a few
seconds each per case. This bounds iteration work, not total elapsed time.
A fit that reaches the cap retains status `ok`, its actual convergence flag,
and its KKT residual. The p ≤ 500 calibration and timed fits retain their
20,000-iteration limit. No large-p calibration sweep is performed.
Calibration results are cached separately. Because workers visit large p
first, each computes a missing p = 500 prerequisite on demand. Workers can
duplicate this deterministic prerequisite calculation; distinct temporary
files and atomic replacement keep the shared cache complete.

The CliqueTrees timed call constructs
`A = sparse(Symmetric(S .* pattern))`, where `pattern` is E with its diagonal
set to true (the dimensioned form of `E .| I`), constructs `ChordalCholesky(A)`, calls
`complete!(F, Symmetric(A))` with its defaults, forms
`F.P' * (sparse(F.L) * sparse(F.L)') * F.P`, and projects its precision
matrix onto E. Its largest absolute off-graph entry is measured before that
projection. Iterations and convergence are unavailable (`NaN`, `missing`).
Exceptions, including `PosDefException`, remain visible.

Before timing CliqueTrees at p ≥ 1000, construct A and
`F = ChordalCholesky(A)` without numerical factorization. Count the triangular
embedding entries as `(ndz(F) + p) ÷ 2 + nlz(F)`, using CliqueTrees' block-size
accessors. This matches its `sparse(F.L)` conversion; `nnz(F.L)` includes the
unused halves of dense diagonal blocks and would overcount. The fill fraction
is this count divided by `p(p+1)/2`. When the fraction exceeds 0.5, skip
CliqueTrees: Newton's method on a nearly dense embedding would defeat the
runtime budget. Record one row with status
`skipped: dense chordal extension (fill <fraction>)`, NaN measurements, and
missing convergence. Skips are separate from failures and excluded from
hypothesis comparisons. Every CliqueTrees row records `embedding_fill`,
including a skipped row or a later exception; it is NaN when not computed.
Other methods also have NaN in that column. A symbolic-construction exception
is recorded as a failure with untimed measurements.

`REFIT_FILL_SKIP_P` sets the smallest dimension subject to the fill rule
(default 1000). Its value is saved in metadata. Setting it to 30 exercises
the rule on the quick cases without changing the graph or the 0.5 threshold.

Every repetition records timing, allocated bytes, settings, iterations,
convergence, KKT residual, objective, and exception type plus its first message
line. `RefitInfo.kkt` supplies iterative package certificates. ChordalCompletion
and CliqueTrees receive an additional certificate outside the timed call:
the largest `abs(inv(Ω)[i,j] - S[i,j]) / sqrt(S[i,i]*S[j,j])` on the graph and
diagonal, computed by Cholesky. The objective is `tr(SΩ) - logdet(Ω)` and is
also evaluated outside timing. A nonpositive definite estimate has objective
`NaN` and KKT residual `Inf`. Package estimates have exact zeros off the graph.

## Timing and restart behavior

Use one Julia thread and one BLAS thread. Each process warms up all five methods
on a p = 8 band problem before any timed case. Each method first runs once.
Its first time determines its total repetitions: 30 below 0.1 seconds, 10
below 1 second, at least 3 with a 30-second timing budget below 60 seconds,
and one at or above 60 seconds. Large-p ADMM uses one fit as specified above.
The medium-duration count is initially
`max(3, floor(30/first_time))`; it stops sooner after three repetitions when
the elapsed fit time plus the estimated next fit exceeds 30 seconds.
The minimum of three takes precedence over that budget. A failed method stops
after its exception. Repetition rounds rotate their starting method through
the fixed five-method order. `GC.gc()` precedes every timed fit. `@timed`
records elapsed seconds and allocated bytes; the garbage collection before
the call is excluded. Certificate checks add untimed work.

Each complete case is atomically written to `cases/<case-id>.tsv` through a
worker-specific `*.tmp` file and skipped on restart. A partial case runs again;
temporary files are ignored. QUIC supports are atomically cached under
`graphs/`, and ADMM calibration under `calibration/`. Workers traverse the same
list in descending p, then ascending n and graph order. Worker k/K takes
entries k, k+K, k+2K, and so on. Keep K fixed when restarting concurrent workers.
Use a fresh output directory when changing code, package versions, or settings;
existing cases and caches are trusted without validation.

`meta.tsv` (or `meta-k.tsv` with `--worker`) records the UTC start date, Julia
and dependency versions, package git commit when available, CPU model class,
BLAS vendor/configuration, and thread counts. Existing metadata is retained
on restart. No machine identity is recorded.

## Run

From the repository root, instantiate the environment using the workspace
dependency helper:

```sh
julia scripts/instantiate.jl benchmark/refit
julia --threads=1 --project=benchmark/refit benchmark/refit/run.jl results/refit
julia --project=benchmark/refit benchmark/refit/summarize.jl results/refit
```

For independent workers sharing an output directory, launch each command in
a separate process, for example:

```sh
julia --threads=1 --project=benchmark/refit benchmark/refit/run.jl results/refit --worker 1/2
julia --threads=1 --project=benchmark/refit benchmark/refit/run.jl results/refit --worker 2/2
```

After all workers finish, run the summarizer once. It prints and writes
`summary.md`, with per-p median/minimum timing tables, objective gates, strict
timing hypotheses H1–H3, failures, skips, nonconvergence, unmatched calibrations, and
the largest CliqueTrees off-graph entry. The pairwise gate checks objective
agreement within 1e-8 among successful fits with KKT ≤ 1e-6. A second gate
rejects any successful, explicitly converged return whose objective exceeds
the best finite successful objective on the case by more than 1e-8 (or is
nonfinite). Cases with fewer than two qualifying methods are listed; their
pairwise difference is defined as zero. Missing hypothesis comparisons are
reported as not evaluated completely. Timing cells show `skipped (fill x)`
for skipped methods and exclude their NaN times. Skipped opponents do not
create unavailable hypothesis comparisons.

Any successful repetition with `converged = false` makes its method's median
time a lower bound for hypothesis evaluation. Write a for the subject's
median and b for an opponent's. With only the opponent capped, `a < b`
supports the hypothesis: the subject beats even the opponent's lower bound.
When `a ≥ b`, the comparison is `unavailable (lower bound)`; a capped
opponent cannot produce a counterexample. With only the subject capped,
`a ≥ b` is a counterexample because the subject is no faster even at its
lower bound; `a < b` is unavailable (lower bound). If both are capped, the
comparison is unavailable (lower bound). With neither capped, `b ≤ a`
is a counterexample and `a < b` supports the hypothesis. H3 applies these
rules to every opponent; a counterexample takes precedence over other
unavailable comparisons. Each hypothesis
lists and counts cases unavailable because of lower bounds, separately from
other unavailable comparisons.

Gates and hypotheses cover completed
cases only, so inspect completeness before drawing full-suite conclusions.

For a small check:

```sh
julia --threads=1 --project=benchmark/refit benchmark/refit/run.jl results/refit-quick --quick
julia --project=benchmark/refit benchmark/refit/summarize.jl results/refit-quick
REFIT_FILL_SKIP_P=30 julia --threads=1 --project=benchmark/refit benchmark/refit/run.jl results/refit-quick-fill --quick
julia --project=benchmark/refit benchmark/refit/summarize.jl results/refit-quick-fill
```

Quick mode uses p = 30, 60, n = 2p, and band, grid (5×6 and 6×10), and ER
graphs. It uses one repetition per applicable method, ADMM tolerance 1e-8
with 1000 iterations, and no calibration. H2 has no eligible quick cases.
Quick timings validate the harness and do not establish full-suite performance.

## Results

Results will be added from the archived run.
