# Random probes and iterative index refinement

## What average probe quality means

Let T=beta'Z have unit variance, with Z standard Gaussian and conditionally
mean-zero X. For a probe q(X), measure its standardized first-order signal by

$$R(q)=\frac{\operatorname{Cov}\{T,q(X)\}^2}{\operatorname{Var}\{q(X)\}}.$$

This is the squared correlation between the probe and the true index. It is
scale invariant. Let D=E[Sigma'(T)]. The Stein loading for a coordinate square
is D[j,j], for a coordinate product is D[j,k], and for two independent random
unit directions a,b it is a'Db. Directions are drawn uniformly on the sphere
and are independent of data and of each other. No probe is selected using D.

For fixed a,b the probe variance is

$$\operatorname{Var}\{(a^TX)(b^TX)\}
=E[(a^T\Sigma(T)a)(b^T\Sigma(T)b)+2(a^T\Sigma(T)b)^2]
-(a^TE[\Sigma(T)]b)^2.$$

This includes both conditional Gaussian noise and variation of the conditional
probe mean. Comparing raw derivative size without this variance can mislead.

With independent spherical directions,

$$E_{a,b}(a^TDb)^2=\frac{\|D\|_F^2}{p^2}.$$

A uniform coordinate square has average squared loading sum_j D[j,j]^2/p.
A uniform distinct coordinate pair has average squared loading
(sum_{j,k} D[j,k]^2-sum_j D[j,j]^2)/(p(p-1)). Neither dominates the dense probe
for every D. In particular, diagonal variation can be diluted by dense probes.

R/(1-R) is also reported. In this Gaussian design it is the signal-to-residual
variance ratio for the linear regression of the probe on Z. It governs the
single-probe unpenalized OLS direction's first-order tangent variance when the
loading is nonzero. Averages of these quantities do not establish average
finite-sample index accuracy. Near-zero loadings and dependence among multiple
probes matter; bank-level information is not a sum of marginal R values.

## Discovery population comparison

Use strongest-amplitude versions of all four approved families, with p=20 and
40. Enumerate every coordinate square and every distinct coordinate pair.
Average 20000 independent random direction pairs and 20000 random squared
projections per configuration. Gaussian quadrature integrates over T. No index
estimation or new precision inference is performed.

At p=40, the mean R values below are multiplied by 1000; larger is better.

| Design | Coordinate square | Coordinate pair | Two random directions | One random direction squared |
| --- | ---: | ---: | ---: | ---: |
| Opposed scale | 51.835 | 0.165 | 3.132 | 3.053 |
| Noncommuting | 0.024 | 3.767 | 4.584 | 4.336 |
| Correlation | 0 | 2.257 | 2.675 | 2.609 |
| Fixed precision diagonal | 1.377 | 2.149 | 2.245 | 3.290 |

The two-direction means have Monte Carlo standard errors 0.030, 0.043, 0.026,
and 0.022 in these table units. Coordinate averages are enumerated, not Monte
Carlo estimates. All values retain the stated quadrature approximation.

Two random directions spread signal more evenly than coordinate products.
In the correlation and opposed-scale families, 97.44% of distinct coordinate
pairs have zero first-order signal. Dense direction pairs almost surely avoid
an exact zero, although their signal can be arbitrarily weak. Their mean R
exceeds a random coordinate pair's mean by approximately 19% in correlation,
22% in noncommuting, and a factor of 19 in opposed scale. However, a coordinate
square has about 16.5 times the dense pair's mean R in opposed scale. In the
fixed-diagonal design at p=20 the dense-pair mean is slightly lower than the
coordinate-pair mean. There is no uniform dominance.

A random squared projection is often comparable to two random directions and
is better in the fixed-diagonal example. A mixed dictionary remains reasonable.
The complete quadratic dictionary spans all these quadratic probes; random
projection banks change compression and weighting, not population information.

Job 13447390 completed in 24 seconds with 49 checks: variance positivity,
correlation bounds, the two spherical loading identities, the Sigma=I variance
identity, and GH121/GH201 quadrature agreement. Compact results and the source
recipe are in `experiments/single_index/results/random_probe_population/`.
Recreate under Discovery Slurm with `random_probe_population.sbatch`, using
the existing coverage-run Manifest and shared Julia 1.13.2 depot. The random
seed is 65193717. No simulation arrays require retention.

## Iterative estimator: original proposal

Treat several quadratic probes as a vector response Y. Its conditional mean
is a vector of one-dimensional curves g(beta'Z). For a coordinate product,
that curve is the corresponding covariance entry. The same construction works
for squared and two-direction projections.

A projected multivariate MAVE refinement would proceed as follows.

1. Initialize beta using the frozen Stein/RRR estimator. Standardize the probe
   responses using the index-training sample. Keep a fixed probe dictionary.
2. At the current projected index, estimate each curve and its slope by local
   polynomial regression. All curves share the same smoothing and derivative
   matrices, so the matrices can be applied to the whole response bank.
3. Freeze the current neighborhood weights and local coefficients. Update beta
   by the common weighted local-linear regression criterion, constrained to unit
   norm. At dimension 20, update all coordinates without a sparsity penalty so
   variables omitted by initial screening can reenter.
4. Recompute projected neighborhoods and repeat. Compare one update with full
   iteration; use held-out or leave-one-out prediction to detect self-fitting.
   Preserve starts, nonconvergence, and adverse outcomes.

One criterion is

$$\sum_{i,j} w_{ij}\left\|Y_j-a_i-s_i b^T(Z_j-Z_i)\right\|_W^2,
\qquad \|b\|_2=1,$$

where a_i and s_i are vectors of local intercepts and slopes. Weights w_ij are
computed from the previous iterate's one-dimensional projected distances and
held fixed within the update. Start with standardized responses and W=I.
Correlated-probe weighting is a separate sensitivity, not an assumed improvement.

This uses changing local slopes rather than only their global averages. It can
therefore exploit variation that cancels in the first-order Stein moment.
Initialization still matters. Alternation does not resolve the sign ambiguity:
b and -b describe the same family after reversing the link functions.

The closest multivariate reference is Hong-Fan Zhang (2021),
[weighted MAVE with group Lasso for multivariate responses](https://www.sciencedirect.com/science/article/pii/S0047259X21000312).
It combines response curves and accounts for correlated responses. Our quadratic
response errors are heteroscedastic and dependent, so its statistical guarantees
must not be imported without checking their assumptions.
[Pautrel and Portier (2026)](https://proceedings.mlr.press/v306/pautrel26a.html)
provide a computational route based on projected neighborhoods and stochastic
Riemannian optimization. Their method would need adaptation to our vector of
quadratic responses. Classical projected MAVE is the simpler first prototype.

The current KernelSmoothing.jl checkout provides `local_polynomial` and
`smoother_matrix`, both with derivative options. A single smoother matrix and
a single slope matrix can be reused across probes. Their construction avoids
smoothing in the full 20-dimensional covariate space.

The first proposed comparison is Stein alone versus one refinement versus
iterated multivariate MAVE, using identical observations and probe banks.
Evaluate direction error before revisiting precision inference. Keep the
precision-estimation sample independent, as in the completed coverage study.
The implementation and newly authorized Discovery comparison are documented in
[MAVE_REFINEMENT.md](MAVE_REFINEMENT.md). This section records the original
proposal; the implementation uses an unconstrained update followed by
normalization and explicitly rescales the local slopes.
