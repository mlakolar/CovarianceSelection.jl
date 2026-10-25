# Fixed-rule aggregation of Gaussian Stein probes

## Purpose and authorization

The user authorized this comparison and Discovery execution on 2026-10-25,
including the compute needed for timely completion. No new precision-coverage
experiment is included. These are experimental estimators, without a new
inference guarantee. All numerical fitting and reporting run under Slurm.

The previous comparison already aggregates probes. The present study asks whether
accounting for probe correlation and sample covariate geometry improves that
aggregation. It also checks whether the design itself is too difficult.

## Frozen design

The index dimension is 20. Its four nonzero coordinates are 2:5; coordinate 1 is
zero. Both Gaussian fixed-support families from MULTIPROBE_INDEX.md are reused.
The chain precision matrix and independent-pairs precision matrix retain the
same support for all index values. X is conditionally Gaussian, and Z is standard
Gaussian. Main index coefficients are equal after normalization.

There are 28 scenarios:

- Each design at response dimensions 20, 40, 160 and n=200, 800, 3200 (18).
- Each design at response dimension 40 and n=400, 1600 (4).
- Each design at response dimension 40, n=800, and signal multipliers 0.5, 1.5 (4).
- Each design at response dimension 40, n=800, with mixed-sign index coefficients (2).

A six-dataset canary includes the largest response dimension and the largest
complete dictionary. Development uses 100 replicates per scenario (2,800 datasets).
Independent confirmation uses 500 (14,000 datasets). Plans and all fixed rules
are generated before observing either stage. Every rule advances to confirmation;
there is no truth-based selection, exact-support gate, or CV search.

Each dataset supplies the same probes to all its methods: 60 sparse projection
squares; mixed probes containing all coordinate squares, 30 random products and
30 projection squares; and the complete quadratic dictionary at response dimensions
20 and 40. Full covariance weighting of 12,880 probes at response dimension 160
is deliberately excluded. That setting compares the two compact dictionaries.
Probe seeds are independent of the observations. Development, confirmation, and
canary use disjoint seed namespaces, recorded in their TOML plans.

## Estimators

For standardized probe responses V and centered covariates Z, define
M=Z'V/(n-1), G=Z'Z/(n-1), and C=V'V/(n-1). The population first moment is beta a'.
Write G=LL'. All directions below are normalized to unit length.

1. Current moment aggregation: leading eigenvector of MM'.
2. Current multivariate OLS aggregation: leading eigenvector of BB', B=G^{-1}M.
3. Rank-one reduced-rank regression (RRR): take the leading eigenvector v of
   L^{-1}MM'L^{-T}, then return L^{-T}v. This is the generalized eigenproblem
   MM'b=lambda G b. It differs from applying SVD to the OLS slope matrix.
4. Covariance-weighted OLS and RRR: replace the identity between the two moment
   matrices by W=((1-rho)C+rho I)^{-1}, for fixed rho=0.05, 0.2, 0.5.
   This downweights redundant probe combinations. It is a feasible candidate,
   not a claim of optimal GMM weighting. Probe covariance need not equal the
   full Stein-moment noise covariance.
5. Learned combined probe: split the observations into two halves. On the training
   half, estimate an RRR direction u and loading a=M'u. Set w=((1-0.2)C+0.2I)^{-1}a,
   normalized to empirical probe variance one. Evaluate the scalar response Vw
   and its first Stein moment on the other half, using training centers/scales.
   Swap halves. Align the two training directions by their inner product and
   average the correspondingly oriented held-out moment vectors. Normalize once
   at the end. No true direction is used for construction or orientation.

Methods 1--4 use either all 20 index variables or joint OLS row-norm screening
that retains 10 variables. OLS and RRR are refit on those variables; direct
moment aggregation restricts the moment vector to those variables. The selected
row set and standardized response dictionary are held constant across weighting
rules. Learned-probe aggregation uses all 20 variables, avoiding confounding by
an additional selection rule. All procedures are non-iterative and do not fit a
nonparametric conditional-mean function.

Two explicitly infeasible diagnostic controls use q(X)=X'[-Omega'(0)]X, an
informative probe constructed from the known simulation precision function.
They fit scalar-response OLS on either all 20 index variables or the four true
variables. They assess intrinsic difficulty and the cost of estimating support;
they are never used for feasible tuning or selection.

## Outcomes and validation

Primary loss is 1-(beta_hat'beta)^2. Report its mean, Monte Carlo standard error,
median, and 90th percentile. Report sign-aligned prediction error, displacement,
and paired loss differences against the matching unweighted method. Failed fits
remain in comparisons with direction loss one and a separate failure record.

As cheap downstream diagnostics, compare displacement with h=n^{-1/5}. On 1,000
independent standard Gaussian evaluation covariates, compare normalized Gaussian
kernel weights at index zero using the estimated and true direction. Report their
L1 distance at bandwidths 0.75h, h, and 1.5h. These are neighborhood diagnostics,
not precision bias, coverage estimates, or a universal accuracy cutoff.

Measure fitting time and separate shared covariance/eigendecomposition setup
costs. Shared setup must not be omitted from computational conclusions.
Validate noiseless identification, generalized-eigen equivalence, orthogonal
probe invariance for unweighted RRR, exact screened zeros, finite learned-probe
outputs, full expected record identities, paired counts, and retained failures.

## Execution and retention

Run identity: 20261025-sig-aggregation. Existing Julia 1.13.2 environment and
Manifest are reused. The canary requests two CPUs / 4 GB. Production requests
64 single-threaded CPUs across four nodes / 128 GB, with a two-hour limit per
stage and GNU Parallel launching exclusive srun steps. Development and confirmation
run sequentially. Use observed canary throughput to adjust resource layout or
wall time without changing frozen scientific settings. Resume only missing atomic
shards; preserve completed failures and original seeds.

Keep the source revision/snapshot, plans, raw compact CSV shards on Discovery,
and version compact summaries, paired results, receipts, and the execution recipe.
Do not advance to precision coverage automatically. Evaluate aggregation gains and
the informative-probe controls first.

## Results and execution status

All stages completed on Discovery. Canary `13432568` passed 16 algebra checks
and 332 expected records. Development `13433102` validated 2,800 datasets and
161,600 records. Independent confirmation `13433104` validated 14,000 datasets
and 808,000 records. There were no duplicate identities or failed fits. Record
counts include shared-setup timing rows. Development and confirmation fitting-input
hash files are identical. First reporting job: `13433770`; figure-layout revision:
`13434797`. No fits were repeated for figure changes.

Representative fresh-confirmation results at n=3,200 and response dimension 40
are below. Entries are mean squared-sine error over 500 paired repetitions.
The weighting level rho=0.2 is one prespecified comparison, not a fitted optimum.

| Design and probes | OLS, all 20 | RRR, all 20 | Weighted RRR, all 20 | Learned combined probe |
|---|---:|---:|---:|---:|
| Chain, mixed | 0.2550 | 0.2520 | 0.3247 | 0.4414 |
| Pairs, complete | 0.1883 | 0.1563 | 0.0801 | 0.2927 |

Relative to unweighted RRR, the paired weighted-RRR difference is +0.07266
(Monte Carlo SE 0.00457) for chain/mixed and -0.07627 (SE 0.00220) for
pairs/complete. Thus accounting for correlation helps in one case and hurts
in the other. For pairs/complete, all-coordinate weighted-RRR loss at fixed
rho=0.05, 0.2, 0.5 is respectively 0.0736, 0.0801, 0.0990. For chain/mixed
it is 0.3461, 0.3247, 0.2896. There is no universally preferred weighting level.

Screening/refitting remains useful. At rho=0.2, pairs/complete RRR retaining
10 variables improves from 0.1224 unweighted to 0.0783 weighted. Chain/mixed
unweighted RRR retaining 10 has loss 0.1997; weighting worsens it to 0.2432.
This particular learned-probe construction does not establish an advantage
over the better direct aggregations. It is not evidence against all learned
or cross-fitted probe constructions.

### Diagnostic controls and downstream neighborhoods

At n=3,200 and response dimension 40, informative-probe OLS on all 20 index
variables has loss 0.0819 for chain and 0.0239 for pairs. Restricting that OLS
fit to the true four variables gives 0.0134 and 0.0037. These probes use the
known simulation precision derivative and are explicitly infeasible controls.
They show that the data permit substantially more accurate index estimation
than the feasible methods considered here. They do not establish a minimax
benchmark or a general guarantee at smaller sample sizes.

At response dimension 160 and n=3,200, informative-probe OLS gives losses
0.0225 (chain) and 0.0067 (pairs). Mixed-probe unweighted RRR retaining 10
variables gives 0.0447 for chain, whereas compact-probe methods still perform
poorly for pairs. The complete dictionary with full covariance weighting was
not tested at response dimension 160. That omission limits conclusions about
high-dimensional response designs.

At response dimension 40, the mean displacement/bandwidth ratio for rho=0.2,
all-coordinate weighted RRR is 2.963 for chain/mixed and 1.417 for pairs/complete,
using h=n^{-1/5}. Their mean kernel-weight L1 distances are 1.183 and 0.765.
The informative-probe, true-support controls have ratios 0.538 and 0.283,
with L1 distances 0.324 and 0.175. These quantify neighborhood changes, not
precision bias or confidence-interval coverage.

Compact tables, all fixed shrinkage levels, paired comparisons, receipts, and
standalone figures are in
[`experiments/single_index/results/aggregation`](../../experiments/single_index/results/aggregation).
The authoritative execution record is
`SingleIndexGraph/.slurm/runs/20261025-sig-aggregation.json`.
