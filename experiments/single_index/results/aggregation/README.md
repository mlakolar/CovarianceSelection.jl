# Fixed-rule Stein probe aggregation

Completed on Discovery, run `20261025-sig-aggregation`. See the
[design and interpretation](../../../../docs/single_index/PROBE_AGGREGATION.md).

- Development: 2,800 datasets, 161,600 records, no failed fits or duplicate identities.
- Independent confirmation: 14,000 datasets, 808,000 records, no failed fits or duplicates.
- Counts include shared setup-timing records. All scientific failures would have been retained.
- Canary: 16 algebra checks, six datasets, all 332 expected records, no failed fits.
- Each summary includes all fixed rules, shrinkage levels, support sizes, and scenarios.
- Paired tables compare each weighted rule with its matching unweighted baseline.
- Figures highlight prespecified rho=0.2; CSVs retain all three levels.
- Neighborhood distances use independent Gaussian evaluation covariates; no precision inference ran.

The full quadratic dictionary is included only for response dimensions 20 and 40.
Response dimension 160 compares compact projection and mixed dictionaries.
Informative-probe and true-support controls use simulation knowledge and are not
feasible procedures. Lower error for these controls does not establish a minimax rate.

## Recreation

Fitting source commit: `0f0c750`. Discovery root:
`/project2/mkolar_1314/SingleIndexGraph/20261025-sig-aggregation`.
The immutable `source.tar.gz`, `plan-source.sha256`, reused Manifest.toml, and
per-stage `inputs-*.sha256` identify the numerical source and plans.
Raw atomic CSV shards remain in that root, outside Git. Compact CSV summaries,
paired tables, receipts, and PDF/PNG figures are versioned here.

Use Julia 1.13.2 and the recorded existing depot. Submit
`sbatch experiments/single_index/aggregation.sbatch canary`, then `development`,
then `confirmation`. Production uses 64 single-threaded workers over four nodes.
Set dependencies as recorded in the run registry. Run
`sbatch experiments/single_index/probe_weighting_report.sbatch` for figures.
The later figure-layout change did not rerun any fitting.

Jobs: canary `13432568`; development `13433102`; confirmation `13433104`;
reporting `13433770`; figure-layout refinement `13434797`. All completed.
