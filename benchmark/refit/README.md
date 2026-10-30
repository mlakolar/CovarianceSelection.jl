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

These are the results of one run on one machine, dated 2026-10-28.
CovarianceSelection 0.1.0 at `d61f21a` used ProximalBase 0.5.4 (master,
`702ba7e`), ProximalSolvers 0.3.0 (master, `f74aeaf`), CliqueTrees 1.19.7,
Krylov 0.10.12, and Julia 1.13.3. The hardware class was an AMD EPYC 9754
128-Core Processor. Four single-threaded workers each used one BLAS thread.
The machine was shared with other workloads. The run covered 40 cases and
2,550 timed repetitions in 30 minutes.

Both the pairwise and converged-objective gates passed. The largest pairwise
objective difference was 5.5e-12. The independent recomputation agreed with
the summary. Objective agreement does not imply equal KKT accuracy: ADMM's
calibration never reached its target. Every ADMM fit used
`abstol = reltol = 1e-10`, with KKT residuals from 2.2e-9 to 7.7e-8, while
DualCoordinate reached ≤ 1e-10, PrimalCoordinate ≤ 6.4e-12, and
ChordalCompletion ≤ 6.5e-16. ADMM's times therefore describe less accurate
fits than the other package methods.

### Median times

Times are in seconds, to three significant digits. Failure and skip labels
are retained from the full summary; — means the method is inapplicable.
The full tables, including minimum times and iterations, are in the archive's
`results/summary.md`.

#### p = 500

| graph, n | ChordalCompletion | DualCoordinate | PrimalCoordinate | RefitADMM | cliquetrees |
|---|---:|---:|---:|---:|---:|
| band, 250 | 0.00252 | 0.0456 | 0.460 | 1.12 | 0.00323 |
| band, 1000 | 0.00256 | 0.0411 | 0.435 | 0.823 | 0.00325 |
| er, 250 | — | 0.0363 | 0.381 | 0.877 | 0.0289 |
| er, 1000 | — | 0.0311 | 0.325 | 0.687 | 0.0192 |
| glasso-dense, 250 | — | 1.51 | 2.81 | 4.99 | fails: PosDefException |
| glasso-dense, 1000 | — | 0.987 | 1.70 | 1.84 | fails: PosDefException |
| glasso-sparse, 250 | — | 0.199 | 0.680 | 7.64 | 1.45 |
| glasso-sparse, 1000 | — | 0.146 | 0.686 | 3.53 | 0.824 |
| grid, 250 | — | 0.0333 | 0.324 | 0.800 | 0.0151 |
| grid, 1000 | — | 0.0320 | 0.380 | 0.759 | 0.0134 |

#### p = 2000

| graph, n | ChordalCompletion | DualCoordinate | PrimalCoordinate | RefitADMM | cliquetrees |
|---|---:|---:|---:|---:|---:|
| band, 1000 | 0.0922 | 0.881 | 49.0 | 45.8 | 0.0589 |
| band, 4000 | 0.0404 | 0.721 | 39.5 | 36.4 | 0.0407 |
| er, 1000 | — | 0.818 | 42.4 | 43.9 | 0.620 |
| er, 4000 | — | 0.641 | 30.4 | 36.7 | 0.394 |
| glasso-dense, 1000 | — | 89.1 | 198 | 178 | skipped (fill 0.940789) |
| glasso-dense, 4000 | — | 60.8 | 117 | 103 | skipped (fill 0.947701) |
| glasso-sparse, 1000 | — | 4.64 | 90.8 | 323 | 133 |
| glasso-sparse, 4000 | — | 3.18 | 55.1 | 144 | 83.3 |
| grid, 1000 | — | 0.799 | 42.5 | 42.6 | 0.122 |
| grid, 4000 | — | 0.643 | 32.4 | 39.7 | 0.0719 |

### Hypotheses and method choice

- **H1: refuted.** PrimalCoordinate was not faster than RefitADMM on every
  case with maximum degree ≤ 50. The counterexamples were the band cases
  at p = 2000, n = 1000 and n = 4000. ADMM attained lower accuracy, so these
  are not comparisons at matched KKT residuals.
- **H2: confirmed on 14 cases.** Wherever DualCoordinate succeeded, it was
  faster than PrimalCoordinate at maximum degree ≤ 10 and p ≥ 500.
- **H3: refuted by one near tie.** On the band case p = 2000, n = 1000,
  ChordalCompletion took 30 repetitions with minimum 0.051 s, median 0.092 s,
  and maximum 0.095 s. CliqueTrees had minimum 0.046 s and median 0.059 s.
  Shared load and this spread limit what can be inferred from the near tie.

DualCoordinate was the fastest iterative method in every case. It never
broke down: no node had at least n neighbours. Thus this run did not exercise
its breakdown regime. PrimalCoordinate remains the fallback when it does
break down. On the band, ER, and grid cases at p = 2000, a primal sweep cost
5.1–6.1 s, making its full fits about as slow as ADMM. RefitADMM provides an
independent check, with achieved accuracy reported alongside its timing.

CliqueTrees was fastest on grids at p ≥ 500, running 2.2–9.0× faster than
DualCoordinate. It was about even with DualCoordinate on Erdős–Rényi graphs,
but about 6–29× slower on the sparse graphical-lasso graphs at p ≥ 500.
It raised `PosDefException` on all four dense graphical-lasso cases at
p ≤ 500, where the package methods agreed on the estimate. Those failures
are failures of this reference call, not evidence that the estimate does
not exist. At p ≥ 1000 the dense graphical-lasso cases were skipped because
their embedding fill fractions were 0.91–0.95.

### Reproduce the archived summary

The [archived run](https://github.com/mlakolar/CovarianceSelection.jl/releases/download/benchmark-evidence/refit-benchmark-20261028-d61f21a.tar.gz)
contains the per-repetition rows, scripts, lockfile, and independent check.
Download it, verify its SHA-256, and extract it:

```sh
curl -fL -o refit-benchmark-20261028-d61f21a.tar.gz \
  https://github.com/mlakolar/CovarianceSelection.jl/releases/download/benchmark-evidence/refit-benchmark-20261028-d61f21a.tar.gz
printf '%s  %s\n' \
  '07f875b17c9557232f00a3c3e6ba10067d24e40fde1072990c761bcc4766b038' \
  'refit-benchmark-20261028-d61f21a.tar.gz' | shasum -a 256 -c -
tar -xzf refit-benchmark-20261028-d61f21a.tar.gz
cd refit-benchmark-20261028-d61f21a
julia --project=scripts scripts/summarize.jl results
```

Alternatively, recompute the checks with the archived independent script:

```sh
python3 -I check/check_benchmark.py results
```
