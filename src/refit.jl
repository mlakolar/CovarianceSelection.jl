"""
    RefitMethod

Supertype of the methods accepted by [`covsel_refit`](@ref):
[`ChordalCompletion`](@ref), [`DualCoordinate`](@ref),
[`PrimalCoordinate`](@ref) and [`RefitADMM`](@ref).
"""
abstract type RefitMethod end

"""
    DualCoordinate(; tol=1e-10, maxiter=1_000)

Maximum likelihood on a known graph by coordinate updates of the covariance
matrix. Each node with neighbours `N` costs O(|N|³ + p|N|) operations per
sweep. The method cannot break down when `S` is positive definite, but can
break down when `S` is singular; in that case it throws an error suggesting
[`PrimalCoordinate`](@ref), even if the maximum likelihood estimate exists.

- `tol`: positive, finite tolerance of the largest change in correlation units
  during a sweep, with an effective minimum of `100eps(T)` for element type `T`.
  Convergence also requires the scaled covariance residual to be at most
  `100tol`.
- `maxiter`: maximum number of sweeps, at least one. Failure to converge gives
  a warning and returns the last estimate.
"""
Base.@kwdef(
  struct DualCoordinate <: RefitMethod
    tol::Float64 = 1e-10
    maxiter::Int = 1_000

    function DualCoordinate(tol::Real, maxiter::Integer)
      tol = Float64(tol)
      isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be positive and finite"))
      maxiter >= 1 || throw(ArgumentError("maxiter must be at least 1"))
      new(tol, maxiter)
    end
  end
)

"""
    PrimalCoordinate(; tol=1e-10, maxiter=1_000)

Maximum likelihood on a known graph by exact minimization over each column
of the precision matrix, keeping its inverse by two rank-one updates. Each
column subproblem is strictly convex and bounded below, so the objective is
non-increasing and the method never breaks down in exact arithmetic, including
when `S` is singular. A numerical Cholesky failure raises `PosDefException`.
It keeps two dense p × p iterates and costs O(p² + |N|³) operations per node
with neighbours `N`, per sweep.

- `tol`: positive, finite tolerance of the largest change in the precision
  matrix in correlation units during a sweep, with an effective minimum of
  `100eps(T)` for element type `T`. Convergence also requires the scaled
  covariance residual to be at most `100tol`.
- `maxiter`: maximum number of sweeps, at least one. Failure to converge gives
  a warning and returns the last estimate. When the maximum likelihood estimate
  does not exist, the objective decreases without bound.
"""
Base.@kwdef(
  struct PrimalCoordinate <: RefitMethod
    tol::Float64 = 1e-10
    maxiter::Int = 1_000

    function PrimalCoordinate(tol::Real, maxiter::Integer)
      tol = Float64(tol)
      isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be positive and finite"))
      maxiter >= 1 || throw(ArgumentError("maxiter must be at least 1"))
      new(tol, maxiter)
    end
  end
)

"""
    RefitADMM(options=ADMMOptions())

Maximum likelihood on a known graph by ProximalSolvers' `ADMMSolver` on
`LogDetLoss(S)` and `IndPattern(E)`, leaving the diagonal free. An iteration
costs one eigendecomposition. `options` controls the residual stopping tests,
iteration limit and penalty parameter, with the data scaling of
[`ADMMOptions`](@ref). The returned estimate is the projected iterate, with
exact zeros outside the graph. Convergence requires the ADMM residual tests
and a positive definite estimate; otherwise it warns and returns the estimate.
"""
struct RefitADMM <: RefitMethod
  options::ADMMOptions
end

RefitADMM() = RefitADMM(ADMMOptions())

"""
    RefitInfo()

Filled in by [`covsel_refit`](@ref) when passed as the `info` keyword:
`iterations`, whether the method `converged`, the largest scaled `change`
during the last sweep, and `kkt`, the largest covariance residual on the graph
and diagonal, in correlation units. `kkt` is `Inf` if the returned precision
matrix is not positive definite. ADMM records `change = NaN`.

The closed form [`ChordalCompletion`](@ref) records `iterations = 0`,
`converged = true`, and `change = kkt = NaN`, without computing an inverse.
"""
mutable struct RefitInfo
  iterations::Int
  converged::Bool
  change::Float64
  kkt::Float64
end

RefitInfo() = RefitInfo(0, false, NaN, NaN)

"""
    covsel_refit(S, E, method::RefitMethod; info=RefitInfo())

Refit a Gaussian precision matrix on a known graph by maximum likelihood:

    minimize  tr(SΩ) - log det Ω,  Ω ≻ 0,
    subject to Ω[i, j] = 0 outside E, i ≠ j.

Its inverse agrees with `S` on the graph and diagonal, the maximum-determinant
completion. Supply a method explicitly:

- [`ChordalCompletion`](@ref): a closed form for chordal graphs.
- [`DualCoordinate`](@ref): coordinate updates of the covariance matrix.
- [`PrimalCoordinate`](@ref): exact column minimization of the precision matrix.
- [`RefitADMM`](@ref): ADMM with projection onto the graph.

`S` must be a finite, symmetric floating-point matrix with positive diagonal;
it may be singular. A `Symmetric` wrapper uses its stored triangle. `E` accepts
the forms of ProximalBase's `IndPattern`: a Boolean matrix, a dense or sparse
estimate whose nonzero entries are the edges, a `SymmetricSparseIterate`, or an
`IndPattern` itself. Stored sparse zeros are ignored. The pattern must be
symmetric and the size of `S`; its diagonal is ignored. Neither input is
modified.

Returns a `SymmetricSparseIterate` of the element type of `S`, with exact
zeros outside the graph; use `Matrix(Ω̂)` for a dense matrix. Pass a
[`RefitInfo`](@ref) as `info` for convergence diagnostics.

`ChordalCompletion` throws an `ArgumentError` for a non-chordal graph or a
clique on which `S` is not positive definite. `DualCoordinate` throws an
`ErrorException` on a failed neighbourhood Cholesky factorization or a
non-positive residual variance, suggesting `PrimalCoordinate`. A numerical
Cholesky failure in `PrimalCoordinate` propagates as `PosDefException`.
Iterative methods warn and return their estimate if they do not converge;
the maximum likelihood estimate may not exist.

# Choosing a method

For a chordal graph (test with [`is_chordal`](@ref)), use [`ChordalCompletion`](@ref):
it is exact and needs no iterations. On the benchmark's band graphs it was
about 9–19× faster than `DualCoordinate` and 25–980× faster than `PrimalCoordinate`.
Otherwise try [`DualCoordinate`](@ref), the fastest iterative method in every
benchmarked case, 1.7–56× faster than `PrimalCoordinate` at p ≥ 500. It can
break down with singular `S`, for example when a node has at least n neighbours.
Use [`PrimalCoordinate`](@ref) then: it avoids breakdown in exact arithmetic,
but a sparse-graph sweep costs O(p³); at p = 2000 it was about as slow as ADMM.
Use [`RefitADMM`](@ref) as an independent check. See `benchmark/refit` for
the measurements and accuracy differences from this single-machine run.

# Example

```julia
S = [1.0 0.4 0.1; 0.4 1.0 0.3; 0.1 0.3 1.0]
X, Z, U = zeros(3, 3), zeros(3, 3), zeros(3, 3)
covsel!(X, Z, U, S, 0.2)
info = RefitInfo()
Ω̂ = covsel_refit(S, Z, PrimalCoordinate(); info)  # relaxed graphical lasso
Matrix(Ω̂)
```
"""
function covsel_refit(S::AbstractMatrix{T}, E, method::RefitMethod;
                      info::RefitInfo=RefitInfo()) where {T<:AbstractFloat}
  p = size(S, 1)
  size(S, 2) == p || throw(DimensionMismatch("S must be a square matrix"))
  p > 0 || throw(ArgumentError("S must not be empty"))
  all(isfinite, S) || throw(ArgumentError("S must be finite"))
  issymmetric(S) || throw(ArgumentError("S must be symmetric"))
  all(>(zero(T)), diag(S)) || throw(ArgumentError("the diagonal of S must be positive"))
  P = _refit_pattern(E, p)
  info.iterations, info.converged, info.change, info.kkt = 0, false, NaN, NaN
  SymmetricSparseIterate(_refit(Matrix{T}(S), P, method, info))
end

function _refit_pattern(E, p::Union{Nothing,Int}=nothing)
  P = copy(E isa IndPattern ? E.P : IndPattern(E).P)
  if p === nothing
    ndims(P) == 2 && size(P, 1) == size(P, 2) ||
      throw(DimensionMismatch("E must be a square pattern"))
  else
    size(P) == (p, p) || throw(DimensionMismatch("E must be a p × p pattern, the size of S"))
  end
  for j in axes(P, 1)
    P[j, j] = false
  end
  issymmetric(P) || throw(ArgumentError("the pattern E must be symmetric"))
  P
end

_tol(method::Union{DualCoordinate,PrimalCoordinate}, ::Type{T}) where {T} =
  max(T(method.tol), 100 * eps(T))

function _refit_symmetrize!(Ω)
  for j in axes(Ω, 2), k in j+1:size(Ω, 1)
    Ω[k, j] = Ω[j, k] = (Ω[k, j] + Ω[j, k]) / 2
  end
  Ω
end

function _refit_certificate(Ω::Matrix{T}, S, P) where {T}
  factor = cholesky(Symmetric(Ω); check=false)
  isposdef(factor) || return (Inf, false)
  W = inv(factor)
  r = one(T) ./ sqrt.(diag(S))
  kkt = zero(T)
  for j in axes(P, 2), k in j:size(P, 1)
    if k == j || P[k, j]
      kkt = max(kkt, abs(W[k, j] - S[k, j]) * r[k] * r[j])
    end
  end
  kkt, true
end

function _refit_warn(method, info)
  if !info.converged
    name, k = nameof(typeof(method)), info.iterations
    change, kkt = info.change, info.kkt
    @warn "covsel_refit with $(name) stopped after $(k) iterations without converging; the maximum likelihood estimate may not exist" change kkt
  end
end

function _refit_coordinate_finish!(Ω, S, P, method, info)
  info.kkt, _ = _refit_certificate(Ω, S, P)
  tol = _tol(method, eltype(S))
  info.converged = info.change <= tol && info.kkt <= 100 * tol
  _refit_warn(method, info)
  Ω
end

function _refit_dual_breakdown(j)
  error("DualCoordinate broke down at node $(j): the covariance iterate is not positive definite on its neighbourhood; the maximum likelihood estimate may still exist, try PrimalCoordinate()")
end

function _refit(S::Matrix{T}, P, method::DualCoordinate, info) where {T}
  p = size(S, 1)
  neighbours = [findall(view(P, :, j)) for j in 1:p]
  d = maximum(length, neighbours)
  A, b, w = Matrix{T}(undef, d, d), Vector{T}(undef, d), zeros(T, p)
  W = copy(S)
  r = one(T) ./ sqrt.(diag(S))
  tol = _tol(method, T)
  for iter in 1:method.maxiter
    change = zero(T)
    for j in 1:p
      N = neighbours[j]
      if isempty(N)
        fill!(w, zero(T))
      else
        block, β = view(A, 1:length(N), 1:length(N)), view(b, 1:length(N))
        copyto!(block, view(W, N, N))
        factor = cholesky!(Symmetric(block); check=false)
        isposdef(factor) || _refit_dual_breakdown(j)
        copyto!(β, view(S, N, j))
        ldiv!(factor, β)
        mul!(w, view(W, :, N), β)
      end
      for k in 1:p
        k == j && continue
        change = max(change, abs(w[k] - W[k, j]) * r[k] * r[j])
        W[k, j] = W[j, k] = w[k]
      end
    end
    info.iterations, info.change = iter, change
    change <= tol && break
  end
  Ω = zeros(T, p, p)
  for j in 1:p
    N = neighbours[j]
    β = view(b, 1:length(N))
    if !isempty(N)
      block = view(A, 1:length(N), 1:length(N))
      copyto!(block, view(W, N, N))
      factor = cholesky!(Symmetric(block); check=false)
      isposdef(factor) || _refit_dual_breakdown(j)
      copyto!(β, view(S, N, j))
      ldiv!(factor, β)
    end
    σ² = S[j, j] - dot(view(S, N, j), β)
    σ² > zero(T) || _refit_dual_breakdown(j)
    Ω[j, j] = one(T) / σ²
    for (a, k) in enumerate(N)
      Ω[k, j] = -β[a] / σ²
    end
  end
  _refit_symmetrize!(Ω)
  _refit_coordinate_finish!(Ω, S, P, method, info)
end

function _refit(S::Matrix{T}, P, method::PrimalCoordinate, info) where {T}
  p = size(S, 1)
  neighbours = [findall(view(P, :, j)) for j in 1:p]
  d = maximum(length, neighbours)
  A = Matrix{T}(undef, d, d)
  b = Vector{T}(undef, d)
  w, u = zeros(T, p), zeros(T, p)
  Ω, W = Matrix(Diagonal(one(T) ./ diag(S))), Matrix(Diagonal(diag(S)))
  r = one(T) ./ sqrt.(diag(S))
  tol = _tol(method, T)
  for iter in 1:method.maxiter
    change = zero(T)
    for j in 1:p
      N = neighbours[j]
      copyto!(w, view(W, :, j))
      w[j] = zero(T)
      w₂₂, s = W[j, j], S[j, j]
      diagonal = one(T) / s
      if isempty(N)
        fill!(u, zero(T))
      else
        block = view(A, 1:length(N), 1:length(N))
        x = view(b, 1:length(N))
        for (a, k) in enumerate(N), (c, l) in enumerate(N)
          block[a, c] = W[k, l] - w[k] * w[l] / w₂₂
        end
        factor = cholesky!(Symmetric(block))
        copyto!(x, view(S, N, j))
        ldiv!(factor, x)
        x ./= -s
        diagonal = one(T) / s - dot(x, view(S, N, j)) / s
        mul!(u, view(W, :, N), x)
        c = dot(view(w, N), x) / w₂₂
        for k in 1:p
          u[k] -= w[k] * c
        end
        u[j] = zero(T)
        for (a, k) in enumerate(N)
          change = max(change, abs(x[a] - Ω[k, j]) / (r[k] * r[j]))
          Ω[k, j] = Ω[j, k] = x[a]
        end
      end
      change = max(change, abs(diagonal - Ω[j, j]) / (r[j] * r[j]))
      Ω[j, j] = diagonal
      @inbounds for l in 1:p
        @simd for k in 1:p
          W[k, l] += s * (u[k] * u[l]) - (w[k] * w[l]) / w₂₂
        end
      end
      for k in 1:p
        W[k, j] = W[j, k] = -s * u[k]
      end
      W[j, j] = s
    end
    info.iterations, info.change = iter, change
    change <= tol && break
  end
  _refit_coordinate_finish!(Ω, S, P, method, info)
end

function _refit(S::Matrix{T}, P, method::RefitADMM, info) where {T}
  p = size(S, 1)
  for j in 1:p
    P[j, j] = true
  end
  options = method.options
  s = _data_scale(tr(S) / p)
  solver = ADMMSolver(; ρ=T(options.ρ * s^2), α=T(options.α), adaptive=options.adaptive,
                      abstol=T(options.abstol / s), abstol_dual=T(options.abstol * s), reltol=T(options.reltol))
  Ω = Matrix{T}(I, p, p)
  result = solve!(Ω, LogDetLoss(S), IndPattern(P), solver, ProximalOptions(; maxiter=options.maxiter))
  _refit_symmetrize!(Ω)
  info.iterations, info.change = result.iterations, NaN
  info.kkt, positive = _refit_certificate(Ω, S, P)
  info.converged = result.converged && positive
  _refit_warn(method, info)
  Ω
end
