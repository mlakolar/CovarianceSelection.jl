# Pointwise coverage after feasible index estimation

## Scope

Compare known beta with independently estimated beta using the **basic**
refitted square-root Lasso estimator. No correction is enabled in this primary
comparison. No support gate, penalty CV, bootstrap, or uniform band is used.
Failures count as noncoverage; conditional coverage is also reported.

The approved designs are opposed scale (amplitude 0.6), noncommuting (0.6), and
correlation (0.55), with 40 responses and 20 covariates. Precision sample sizes
are 200, 400, 800, and 1600. Index training initially uses the same sample size,
independently, so estimated-index total sample size is twice the precision n.
The known-index comparator uses the identical precision observations.

## Index estimation

Opposed scale uses the mixed probe bank: all 40 squares, 30 random products,
and 30 squared random projections. The other designs use all 820 quadratic
products including squares. Probe randomness is independent of observations
and truth. Center and standardize probes, regress them on the 20 covariates,
rank covariates by OLS coefficient row norms, retain 10, and refit without
penalty. Aggregate by the rank-one reduced-rank regression direction. An
unscreened 20-covariate sensitivity is reserved for the paired stage.

Orientation uses the fixed external reference vector with entries 1/j. It
never uses true beta or a fixed nonzero first coordinate. Record orientation
errors; never discard them. The primary evaluation point is zero, which is
sign invariant; secondary points are -0.5 and 0.5. The DGP beta has entries
(1,1,-1,-1)/2 on coordinates 2:5 and zero elsewhere.

Training and precision samples have independent response permutations. The
worker aligns training response columns with the precision sample's response
labels before fitting. This preserves the shared observable response identities.
Targets are selected by canonical graph role and mapped through the precision
sample permutation, solely to identify prespecified evaluation entries.

## Precision inference

Use singleton residualization for diagonals and common-complement two-response
residualization for off-diagonals. Positive pilot-kernel weights define weighted
square-root Lasso; refit each selected support by weighted least squares.
Use a degree-two local polynomial on residual products, invert its intercept,
and use HC0 normal intervals. Residual degrees-of-freedom, HC3, and inverse-bias
corrections are all disabled. Singular refits are failures, without ridge repair.

The fixed normal-score penalty uses tail probability 0.1 per regression and
C in {0.75,1.1,1.25}. Its predictor score scale uses squared kernel weights;
no CV is performed. Final bandwidth is hfactor*n^(-1/5), with hfactor in
{1,1.5,2}; pilot bandwidth ratios are {1,2,3}. Units are the unit-variance index.

## Stages and retention

1. Canary: two paired replications at n=200 and 800 in each family, all 27
   tuning combinations, to exercise targets, orientation, and numerical failures.
2. Known-index development: 300 replications, all 12 design/sample-size settings,
   all 27 tuning combinations. This calibrates precision inference independently
   of index error.
3. Freeze a small set using known-index development only. Compare known and
   estimated directions on matched precision samples, then confirm frozen basic
   specifications with 2000 independent replications. Confirmatory seeds are
   disjoint from development. Include unscreened-index and doubled-training
   sensitivities on selected settings.
4. If known-index coverage is poor, investigate corrections in separate arms:
   degrees of freedom only, HC3 only, and inverse-bias only. Keep the basic arm
   as reference. Do not replace the primary procedure silently.

Report each entry and point separately: unconditional/conditional coverage,
availability, interval length, bias, empirical SD, RMS standard error,
studentized error mean/SD, index error, orientation errors, local effective
sample size, and selected nuisance support size. Training uncertainty is not
included automatically in HC0: the experiment tests that approximation.

Numerical workers and summaries run under Slurm on Discovery. Retain bounded
raw CSV shards there, and compact summary/failure tables, receipts, source
snapshot hash, settings, and job identifiers in the run record. No manuscript
or precision-band claim is changed by this experiment.

### Frozen development selection rule

Keep the fixed reference C=1.1, final factor=1.5, pilot ratio=2. Also select the
configuration minimizing mean absolute known-beta coverage error from 0.95,
plus twice mean failure rate, with equal weight on all entry/point/n cells.
If this equals the reference, retain the next-ranked configuration as the
second choice. This is simulation-based development calibration, not a general
data-driven bandwidth selector. Preserve the full grid, including adverse cells.

Trigger the separate known-beta correction experiment if the development-best
specification has coverage below 0.925 or availability below 0.99 in any cell
with n>=400. This is an exploratory trigger, not a multiple-testing claim.
Corrections are tested separately; the 2000-replication confirmation still
uses the frozen basic procedure, regardless of correction results.

The allocation pipeline runs grid, paired development, triggered separate
corrections, basic confirmation, and doubled-training sensitivity. The paired
stage includes the all-20-covariate index fit. Seeds differ by stage; no
estimated-index coverage influences the frozen tuning choices.
