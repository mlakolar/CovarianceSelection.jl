# Index orientation correction

A zero first coefficient is allowed by the single-index model. The previous
default required coordinate 1 to survive thresholding and rejected an otherwise
nonzero direction when that coordinate was absent. This was an implementation
restriction, not an identification requirement of the model.

`estimate_index` and `select_index_tuning` now default to `anchor=nothing`.
Normalize the selected vector to unit length and make its largest-magnitude
component positive, choosing the first maximizer at a tie. A caller may request
an explicit nonzero `anchor=j`, or supply `reference=original_direction` to align
by a positive inner product. The OLS support-refitting helper uses the same rule.
An entirely zero selected vector still fails. No unthresholded or forced-support
fallback is silently substituted for that distinct problem.

Sign is part of the curve's coordinate system: `(direction, curve(t))` and
`(-direction, curve(-t))` describe the same distribution. Bootstrap directions
must be aligned with the original direction. Automatic canonical orientation
is stable at a unique largest population component, or when tied maximizers
have the same sign. Opposite-sign population ties require a stable external
reference or another stated convention. An exactly orthogonal reference cannot
orient the estimate and is rejected.

The historical Gaussian and CV workers request `anchor=1` explicitly. Their
existing seed configurations still represent the old protocol. Earlier normal
coverage and bootstrap availability must be interpreted accordingly. They have
not been recomputed with the corrected convention.

## Bounded paired diagnosis on Discovery

Job 13423220 ran 92 existing inference checks, 107 fixed-support/worker checks,
20 orientation checks and two OLS-refit checks. All passed. The tests cover an
exactly zero first coefficient, reference alignment across opposite-sign near
ties, CV score invariance under sign reversal, explicit-anchor rejection and
empty-support rejection.

The paired diagnosis uses 500 datasets at each of n=800 and n=3200, p=40,d=20.
For each dataset, permuting Gaussian covariate labels moves a truly zero
coefficient to coordinate 1 while preserving the responses, actual index and
all signal magnitudes. Both policies select support with the frozen ordinary-
mean threshold multiplier 1.25 and then refit OLS. The output has 4000 rows,
with unique identities and every failure retained.

| n | True first coefficient | Policy | Available / 500 | Empty support | Missing fixed anchor |
|---:|---:|---|---:|---:|---:|
| 800 | 0.5 | Fixed first | 177 | 78 | 245 |
| 800 | 0.5 | Automatic | 422 | 78 | 0 |
| 800 | 0 | Fixed first | 0 | 78 | 422 |
| 800 | 0 | Automatic | 422 | 78 | 0 |
| 3200 | 0.5 | Fixed first | 499 | 0 | 1 |
| 3200 | 0.5 | Automatic | 500 | 0 | 0 |
| 3200 | 0 | Fixed first | 1 | 0 | 499 |
| 3200 | 0 | Automatic | 500 | 0 | 0 |

There were no other failures. Automatic orientation gives identical availability
and direction-error summaries before and after this relabeling. At n=800, it
does not repair the 78 empty-support outcomes. Its squared direction error,
averaged over successful fits, is about 0.711; at n=3200 it is about 0.018.
These conditional means are not directly comparable with the fixed-first means,
which average over different subsets of datasets. Availability alone does not
establish accurate index estimation or valid precision confidence intervals.

Compact output and the run record are in
`experiments/single_index/results/index_orientation/`. This is a first-stage
diagnosis, not a replacement pointwise or bootstrap coverage study.
