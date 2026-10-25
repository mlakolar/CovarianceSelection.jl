# First-stage index selection

The diagonal/CV confirmation found a remaining generated-index gap even with
oracle nuisance regressions. `select_index_tuning` addresses index estimation
separately. It is an optional empirical extension; it does not replace the
original median-of-means default or establish generated-index coverage.

## Validation objective

The identifying Stein moment has the form `E(V)=c*gamma`, where `V` is the
known covariate score times the prespecified probe response product. The index
direction is unit-normalized with an explicit sign convention. For each
candidate and fold, estimate a direction on the training observations. Score
its inner product with the held-out moment, multiplied by the sign of the inner
product between the training raw moment and training direction. This aligns
validation with the raw training direction, independently of the reported sign
convention. For a fixed positive anchor it reproduces the original score based
on the training anchor moment. For a fixed training fit, the conditional
population objective is
`c*sign(dot(training_moment,direction))*dot(gamma,direction)`.
Under correct raw-moment alignment this is `abs(c)` times the absolute cosine
with the true direction. A weak empirical moment need not align correctly.
The criterion measures direction accuracy up to sign; a separate convention
defines the sign of the reported index axis.

Thus held-out alignment targets direction accuracy without using the true
index, covariance, or empirical confidence coverage. All candidates use the
same folds and validation observations. The selected candidate is refitted
using all observations. The second-stage index is still estimated on the same
data as the precision curve; this is not an honest sample-splitting argument.

## Candidate estimators

`estimate_index(...; unshrink=true)` thresholds the moment to select a support,
then retains the raw moment on that support before normalizing. The default
soft-thresholds its magnitude. Unshrinking can remove directional shrinkage
bias and can amplify retained noise. Both variants remain available to CV.

`select_index_tuning` accepts `block_options`. Five blocks give the original
coordinatewise median-of-means moment; one block gives the ordinary mean.
The default `(blocks,)` preserves the requested moment estimator. The
light-tail refinement compares one and five blocks. This is not a claim that
the ordinary mean has the robustness of median-of-means. The validation mean
itself is sensitive to tails.

The threshold grid is prespecified. Candidate failures receive negative
infinite alignment. Ties prefer smaller support, a larger multiplier,
unshrinking, then fewer blocks. The final full-data fit may still fail; this
failure is retained. A selected probe with a weak or zero Stein constant is
not repaired by this selector.

```julia
selected = select_index_tuning(X, score;
    multipliers=[0.25,0.5,0.75,1.0,1.25], block_options=[1,5],
    unshrink_options=[false,true], rng=MersenneTwister(19))
index = Z*selected.estimate.direction
```

Return values include the final `IndexEstimate`, multiplier, block count,
unshrinking flag, candidate alignments/support sizes/failure messages, and
fold assignment. Simulation CSVs retain selected index multiplier, support
size, unshrinking flag, and block count, in addition to direction error.

The 701000000 namespace develops five-block alignment CV. The separate
719000000 namespace compares one- and five-block candidates. The independent
751000000 confirmation freezes this candidate class before its outcomes; see
[INDEX_CONFIRMATION_PROTOCOL.md](INDEX_CONFIRMATION_PROTOCOL.md) and
[RESULTS.md](RESULTS.md) for the completed assessment. Weak identification, high-dimensional index estimation error,
and heavy-tailed validation remain explicit limitations.
