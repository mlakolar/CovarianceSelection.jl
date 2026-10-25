# Four Gaussian covariance families: index comparison

See the [designs and interpretation](../../../../docs/single_index/FOUR_GAUSSIAN_DESIGNS.md).

- Population diagnostics: 24 configurations, 162 checks.
- Canary: 98 checks, eight datasets, 444 expected records, zero failed fits.
- Comparison: 240 settings, 300 paired repetitions each, 72,000 datasets.
- All 3,996,000 expected records are present, with no duplicate identities.
- Counts include 216,000 shared-setup timing records.
- Fifty-six split-moment fits had no positive empirical signal eigenvalue. All
  receive direction loss one; no scientific failures were discarded.
- Selected-support identities and all 120 d=4 oracle-control matches were checked.
- Canary and production fitting-input hashes agree.
- All numerical work and figures ran on Discovery. No precision estimation or inference ran.

The CSV files retain every amplitude, method, support size, shrinkage level and
sample size. Mean errors, Monte Carlo errors, medians, 90th percentiles, support
diagnostics, and bandwidth-relative neighborhood diagnostics are reported.
Paired contrasts use the same datasets and probes. Figures show prespecified
rho=0.2 comparisons; the CSV files also retain rho=0.05 and 0.5.
Known informative-probe and true-support fits are infeasible diagnostic controls.
They are not inputs to feasible tuning or an information lower bound.

## Recreation

Run: `20261025-sig-four-designs`. Fitting source: commit `008d7ea`.
Discovery root: `/project2/mkolar_1314/SingleIndexGraph/20261025-sig-four-designs`.
Reuse the recorded Julia 1.13.2 Manifest/depot. The complete source and frozen
plans are in the repository; stage input hashes are included here.

First submit `four_design_population.sbatch`. After its checks pass, run
`four_design.sbatch canary`; then run `four_design.sbatch comparison`. Scripts
are under `experiments/single_index/`. Production uses 64 single-threaded workers
over four nodes. Run `four_family_report.sbatch` for figures and
`four_family_review.sbatch` for the targeted failure/support review.

Successful jobs: population 13445376; canary 13445394; comparison 13445513;
reporting 13445520; failure review 13445799; final figure refinement 13445820.
The two population-driver startup failures and their repairs are preserved in the
run registry. They occurred before estimator fitting.

Raw atomic CSV shards remain on Discovery. These files are compact reporting
inputs, receipts, retained failure records and standalone figures.
