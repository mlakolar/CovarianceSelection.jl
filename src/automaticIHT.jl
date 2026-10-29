# IHT uses a single product-space step chosen by backtracking.
const _IHTGradientMethod = Union{ProximalSolvers.ProximalGradientDescent,
  ProximalSolvers.AccProxGradDescent, ProximalSolvers.MonotoneAccProxGradDescent}

function _iht_covariance(Σ::AbstractMatrix)
  p = size(Σ, 1)
  size(Σ, 2) == p || throw(DimensionMismatch("the covariance must be square"))
  p > 0 && issymmetric(Σ) && all(isfinite, Σ) ||
    throw(ArgumentError("the covariance must be nonempty, finite and symmetric"))
  p
end

function _iht_start(::Type{T}, p, s, r, init; covariance = nothing) where {T}
  s >= 0 || throw(ArgumentError("s must be nonnegative"))
  0 <= r <= p || throw(ArgumentError("r must be between zero and p"))
  if init === nothing
    S, U = SymmetricSparseIterate(T, p), zeros(T, p, r)
    for i in 1:p
      if covariance !== nothing
        covariance[i, i] > 0 || throw(ArgumentError("automatic initialization needs positive covariance diagonal entries; supply init=(S,U)"))
        S[i, i] = inv(covariance[i, i])
      end
      i <= r && (U[i, i] = sqrt(covariance === nothing ? one(T) : S[i, i]) / T(10))
    end
  else
    length(init) == 2 || throw(ArgumentError("init must be (S,U)"))
    S0, U0 = init
    size(S0) == (p, p) && size(U0) == (p, r) || throw(DimensionMismatch("init must have shapes (p,p) and (p,r)"))
    issymmetric(S0) || throw(ArgumentError("the initial sparse matrix must be symmetric"))
    S, U = SymmetricSparseIterate(Matrix{T}(S0)), Matrix{T}(U0)
  end
  all(isfinite, S) && all(isfinite, U) || throw(ArgumentError("the initial blocks must be finite"))
  S, U
end

# Dense matrix buffers are reused across value, gradient and line-search calls.
struct _LatentLogdet{T} <: AbstractFunction
  Σ::Matrix{T}
  Θ::Matrix{T}
  factor::Matrix{T}
  inverse::Matrix{T}
end
function _LatentLogdet(Σ::AbstractMatrix{T}) where {T<:AbstractFloat}
  p = _iht_covariance(Σ)
  _LatentLogdet(Matrix(Σ), zeros(T, p, p), zeros(T, p, p), zeros(T, p, p))
end
struct _LatentDTrace{T,F} <: AbstractFunction
  loss::F
  Θ::Matrix{T}
  raw::Matrix{T}
end
function _LatentDTrace(Σx::Symmetric{T}, Σy::Symmetric{T}) where {T<:AbstractFloat}
  p = _iht_covariance(Σx)
  _iht_covariance(Σy) == p || throw(DimensionMismatch())
  _LatentDTrace(direct_difference_loss(Σx, Σy), zeros(T, p, p), zeros(T, p, p))
end
ProximalBase.issmooth(::Type{<:_LatentLogdet}) = true
ProximalBase.issmooth(::Type{<:_LatentDTrace}) = true

function _iht_theta!(Θ, x, σ)
  _densify!(Θ, x[1])
  mul!(Θ, x[2], transpose(x[2]), eltype(Θ)(σ), one(eltype(Θ)))
end
function _iht_factor!(f::_LatentLogdet, x)
  _iht_theta!(f.Θ, x, -1)
  all(isfinite, f.Θ) || return nothing
  F = cholesky!(Symmetric(copyto!(f.factor, f.Θ)); check = false)
  issuccess(F) ? F : nothing
end
function (f::_LatentLogdet{T})(x::Tuple) where {T}
  F = _iht_factor!(f, x)
  F === nothing ? T(Inf) : dot(f.Σ, f.Θ) - logdet(F)
end
function ProximalBase.value_and_gradient!(G::Tuple, f::_LatentLogdet{T}, x::Tuple) where {T}
  F = _iht_factor!(f, x)
  F === nothing && return T(Inf) # no inverse or gradient outside the cone
  fill!(f.inverse, zero(T))
  for i in axes(f.inverse, 1)
    f.inverse[i, i] = one(T)
  end
  ldiv!(F, f.inverse)
  # Triangular solves leave the inverse symmetric only up to rounding; halving first avoids overflow.
  G[1] .= f.Σ .- (f.inverse ./ T(2) .+ transpose(f.inverse) ./ T(2))
  mul!(G[2], G[1], x[2], T(-2), zero(T))
  dot(f.Σ, f.Θ) - logdet(F)
end
(f::_LatentDTrace)(x::Tuple) = f.loss(_iht_theta!(f.Θ, x, 1))
function ProximalBase.value_and_gradient!(G::Tuple, f::_LatentDTrace{T}, x::Tuple) where {T}
  value = value_and_gradient!(f.raw, f.loss, _iht_theta!(f.Θ, x, 1))
  G[1] .= f.raw ./ T(2) .+ transpose(f.raw) ./ T(2)
  mul!(G[2], G[1], x[2], T(2), zero(T))
  value
end

# Collect every accepted outer iterate, so trace display runs at every iteration.
# A legacy callback's return value is ignored; options.callback can stop solve!.
function _automatic_iht!(x, f, g, method, options::ProximalOptions, accepted)
  method isa _IHTGradientMethod || throw(ArgumentError("automatic IHT needs ProximalGradientDescent, AccProxGradDescent or MonotoneAccProxGradDescent"))
  isfinite(g(x)) || throw(ArgumentError("the initial sparse matrix exceeds the sparsity budget"))
  initial = f(x)
  isfinite(initial) || throw(ArgumentError("the initial loss must be finite; for latent log-det use S - U*U' positive definite"))
  fvals = typeof(initial)[initial]
  cb = function (record)
    state = options.store_trace ? last(record) : record
    push!(fvals, state.fg_x)
    accepted === nothing || accepted()
    options.callback === nothing || state.iteration % options.print_every != 0 || return options.callback(record)
    false
  end
  settings = (; (key => getfield(options, key) for key in fieldnames(typeof(options)))...)
  result = solve!(x, f, g, method, ProximalOptions(; settings..., callback = cb, print_every = 1))
  fvals, result.iterations
end

"""
    difference_precision_iht(Σx, Σy, s; method=ProximalGradientDescent(), options=ProximalOptions(), init=nothing)

IHT with an automatically backtracked step and ProximalBase's Frobenius
projection onto at most `s` lower-triangle coordinates, including the diagonal.
`init` optionally supplies a symmetric starting matrix. Returns `(Δ, fvals)`,
with the initial objective and every accepted iteration in `fvals`.
"""
function difference_precision_iht(Σx::Symmetric{T}, Σy::Symmetric{T}, s::Int;
  method = ProximalSolvers.ProximalGradientDescent(), options::ProximalOptions = ProximalOptions(), init = nothing) where {T<:AbstractFloat}
  p = _iht_covariance(Σx)
  _iht_covariance(Σy) == p || throw(DimensionMismatch())
  Δ, _ = _iht_start(T, p, s, 0, init === nothing ? nothing : (init, zeros(T, p, 0)))
  fvals, _ = _automatic_iht!(Δ, direct_difference_loss(Σx, Σy), IndBallL0(s), method, options, nothing)
  Δ, fvals
end

"""
    covsel_latent_iht!(Σ, s, r; method=ProximalGradientDescent(), options=ProximalOptions(), init=nothing, callback=nothing)

Estimate `Θ = S - L`, with `L = UU' ⪰ 0` of rank at most `r`, by iterative
hard thresholding of `S` on `tr(ΣΘ) - logdet(Θ)`. This is the marginal precision
of the observed variables in a Gaussian graphical model with `r` latent
variables (Chandrasekaran, Parrilo and Willsky, 2012). At most `s` off-diagonal
pairs in `S` are nonzero; the diagonal is unconstrained.

Step sizes are chosen by backtracking, and the line search rejects points
outside the positive definite cone. With `MonotoneAccProxGradDescent` the
objective does not increase up to rounding. This is a local method for a
nonconvex problem and need not yield a global minimum.

By default `S = Diagonal(1 ./ diag(Σ))`, `U[i, i] = sqrt(S[i, i]) / 10` for
`i ≤ r`, and the other entries of `U` are zero. The covariance diagonal must
be positive. `init = (S, U)` is copied and must have `S - UU'` positive definite
and respect the sparsity budget.

Returns `(S, L, U, fvals, iter)`, with `S` a `SymmetricSparseIterate`, `L = UU'`,
`fvals` containing the objective at the initial point and every accepted
iteration, and `iter` the number of iterations. `covsel_objective(Σ, S, L)` is
the objective at the result. `callback(Σ, S, L)` observes only accepted iterates;
its return value is ignored. `options.callback` keeps ProximalSolvers' convention.
"""
function covsel_latent_iht!(Σ::AbstractMatrix{T}, s::Int, r::Int;
  method = ProximalSolvers.ProximalGradientDescent(), options::ProximalOptions = ProximalOptions(), init = nothing, callback = nothing) where {T<:AbstractFloat}
  p = _iht_covariance(Σ)
  S, U = _iht_start(T, p, s, r, init; covariance = Σ)
  L = zeros(T, p, p)
  accepted = () -> (mul!(L, U, transpose(U)); callback === nothing || callback(Σ, S, L); nothing)
  fvals, iter = _automatic_iht!((S, U), _LatentLogdet(Σ), SeparableSum(IndBallL0(s; diagonal = false), ZeroFunction()), method, options, accepted)
  mul!(L, U, transpose(U))
  S, L, U, fvals, iter
end

"""
    difference_latent_precision_iht!(Δ, L, U, Σx, Σy, s; method=ProximalGradientDescent(), options=ProximalOptions(), callback=nothing)
    difference_latent_precision_iht(Σx, Σy, s, r; init=nothing, kwargs...)

Automatic-step IHT for the direct difference loss at `Δ+UU'`. The in-place
method starts from and updates `Δ,U`, and overwrites `L=UU'`. Its sparsity budget
counts diagonal coordinates. Returns `(Δ,L,U,fvals,iter)` with initial and every
accepted objective. The return value of `callback(Δ, L)` is ignored; stopping
callbacks belong in `options`. The convenience method copies `init=(Δ,U)`, or
starts from zero `Δ` and a small deterministic orthogonal `U`.
"""
function difference_latent_precision_iht!(Δ::SymmetricSparseIterate{T}, L::StridedMatrix{T}, U::StridedMatrix{T},
  Σx::Symmetric{T}, Σy::Symmetric{T}, s::Int;
  method = ProximalSolvers.ProximalGradientDescent(), options::ProximalOptions = ProximalOptions(), callback = nothing) where {T<:AbstractFloat}
  p = _iht_covariance(Σx)
  _iht_covariance(Σy) == p && size(Δ) == size(L) == (p, p) && size(U, 1) == p || throw(DimensionMismatch())
  s >= 0 && size(U, 2) <= p || throw(ArgumentError("s must be nonnegative and the factor rank at most p"))
  all(isfinite, Δ) && all(isfinite, U) || throw(ArgumentError("the initial blocks must be finite"))
  accepted = () -> (mul!(L, U, transpose(U)); callback === nothing || callback(Δ, L); nothing)
  fvals, iter = _automatic_iht!((Δ, U), _LatentDTrace(Σx, Σy), SeparableSum(IndBallL0(s), ZeroFunction()), method, options, accepted)
  mul!(L, U, transpose(U))
  Δ, L, U, fvals, iter
end
function difference_latent_precision_iht(Σx::Symmetric{T}, Σy::Symmetric{T}, s::Int, r::Int;
  init = nothing, kwargs...) where {T<:AbstractFloat}
  p = _iht_covariance(Σx)
  Δ, U = _iht_start(T, p, s, r, init)
  difference_latent_precision_iht!(Δ, zeros(T, p, p), U, Σx, Σy, s; kwargs...)
end
