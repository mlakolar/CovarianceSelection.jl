# this file implements direct difference estimation for Δ = Ωx - Ωy

####################################
#
# loss function evaluation
#
####################################

"""
    diff_loss(Σx, Δ, Σy, p)
    diff_loss(Σx, Δ, U, Σy, p)

Residual of the estimating equation of direct difference estimation,

    ‖Σx⋅Δ⋅Σy - (Σy - Σx)‖ₚ,

as an entrywise norm: `p = 2` for the Frobenius norm or `p = Inf` for the
largest absolute entry. With the low rank factor `U` (`p × r`) it is the
residual of `Δ + UU'`. `Σx` and `Σy` are `Symmetric` covariance matrices and
`Δ` is a `SymmetricSparseIterate` or a symmetric matrix.

With the population covariances the residual is zero at `Δ = Ωx - Ωy`. With
the sample covariances of held-out data it measures how well an estimate fits
new data, so it serves as a validation loss for choosing the penalty of
[`difference_precision_active_shooting`](@ref).
"""
diff_loss(Σx::Symmetric{T}, Δ::AbstractMatrix{T}, Σy::Symmetric{T}, p::Real) where {T<:AbstractFloat} =
  _diff_loss(Σx * _dense(Δ) * Σy - Σy + Σx, p)

diff_loss(Σx::Symmetric{T}, Δ::AbstractMatrix{T}, U::StridedMatrix{T}, Σy::Symmetric{T}, p::Real) where {T<:AbstractFloat} =
  _diff_loss(Σx * (_dense(Δ) + U * U') * Σy - Σy + Σx, p)

# products with a dense matrix go to BLAS
_dense(Δ::SymmetricSparseIterate) = Matrix(Δ)
_dense(Δ::AbstractMatrix) = Δ

function _diff_loss(R::AbstractMatrix, p::Real)
  if p == 2
    return norm(R)
  elseif p == Inf
    return maximum(abs, R)
  else
    throw(ArgumentError("p should be 2 or Inf"))
  end
end




####################################
#
# loss tr(Σx⋅Δ⋅Σy⋅Δ)/2 - tr(Δ(Σy-Σx))
#
####################################
"""
    direct_difference_loss(Σx::Symmetric, Σy::Symmetric)

The direct difference loss `tr(Σx⋅Δ⋅Σy⋅Δ)/2 - tr(Δ⋅(Σy - Σx))` of a symmetric
`Δ`, as ProximalBase's `KronQuadratic(Σx, Σy, Σx - Σy)`. Its minimizer with
population covariances is `Δ = Ωx - Ωy`, the difference of the precision
matrices (Zhao, Cai and Li, 2014). ProximalSolvers' coordinate descent minimizes
it, plus a penalty, over a `SymmetricSparseIterate`, whose coordinates are the
entries of the lower triangle.
"""
direct_difference_loss(Σx::Symmetric, Σy::Symmetric) = KronQuadratic(Σx, Σy, Σx - Σy)


#####################################################
#
# Direct difference estimation using Active Shooting
#
#####################################################

"""
    difference_precision_active_shooting(Σx, Σy, g, method=CDSolver())
    difference_precision_active_shooting!(x, Σx, Σy, g, method=CDSolver())

Direct estimation of the difference of precision matrices `Δ = Ωx - Ωy`
(Zhao, Cai and Li, 2014) by active shooting coordinate descent:

    minimize  tr(Σx⋅Δ⋅Σy⋅Δ)/2 - tr(Δ⋅(Σy - Σx)) + g(Δ)

where `g` is a separable penalty of ProximalBase on the entries of the
symmetric `Δ`, such as `NormL1(λ)`, or `NormL1(λ, W)` with a symmetric `p × p`
matrix of weights `W`; an off-diagonal entry is counted twice, as `Δ[r, c]`
and `Δ[c, r]`. `method` is ProximalSolvers' `CDSolver`. The in-place version warm
starts from the `SymmetricSparseIterate` `x`. Returns the
`SymmetricSparseIterate` estimate.
"""
function difference_precision_active_shooting!(
  x::SymmetricSparseIterate,
  Σx::Symmetric,
  Σy::Symmetric,
  g::AbstractFunction,
  method::CDSolver=CDSolver())

  solve!(x, direct_difference_loss(Σx, Σy), g, method)
  x
end

difference_precision_active_shooting(
  Σx::Symmetric,
  Σy::Symmetric,
  g::AbstractFunction,
  method::CDSolver=CDSolver()) =
  difference_precision_active_shooting!(SymmetricSparseIterate(float(eltype(Σx)), size(Σx, 1)), Σx, Σy, g, method)

"""
    difference_precision_refit(Σx, Σy, S)

Minimize the direct difference loss `tr(Σx⋅Δ⋅Σy⋅Δ)/2 - tr(Δ⋅(Σy - Σx))` over
symmetric `Δ` supported on `S`, a collection of linear or Cartesian indices into
a `p × p` matrix (an entry and its transpose are the same coordinate). The
minimizer solves the linear system of the Hessian on the support; when it is
singular, a solution of the rank-deficient system.
"""
function difference_precision_refit(Σx::Symmetric, Σy::Symmetric, S)
  p = size(Σx, 1)
  entries = CartesianIndices((p, p))
  support = unique!(sort!([(rc = Tuple(entries[i]); tril_index(p, max(rc...), min(rc...))) for i in S]))
  refit!(SymmetricSparseIterate(float(eltype(Σx)), p), direct_difference_loss(Σx, Σy); support)
end

# D = Δ as a dense symmetric matrix
function _densify!(D::StridedMatrix{T}, Δ::SymmetricSparseIterate{T}) where {T<:AbstractFloat}
  fill!(D, zero(T))
  for k=1:nnz(Δ)
    r, c = tril_cartesian(Δ.p, Δ.nzval2ind[k])
    D[r, c] = D[c, r] = Δ.nzval[k]
  end
  D
end
