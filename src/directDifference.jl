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




########################################################################
#
# Direct difference estimation using iterative hard thresholding
#
########################################################################

# gradient (A + A')/2 + Σx - Σy of the loss with respect to Θ; the IHT solvers
# pass Σx and Σy as dense matrices, which are faster to index than a Symmetric
function diff_precision_grad_delta!(
  grad_out::StridedMatrix{T},
  A::StridedMatrix{T},
  Σx::AbstractMatrix{T}, Σy::AbstractMatrix{T}) where {T<:AbstractFloat}

  p = size(Σx, 1)

  # (a, b) indices of updated grad_out[a,b] element
  for a=1:p, b=a:p
    grad_out[a,b] = (A[a,b] + A[b,a])/2 + Σx[a,b] - Σy[a,b]
    # update the other symmetric element
    if a != b
      grad_out[b,a] = grad_out[a,b]
    end
  end

  grad_out
end


# gradient with respect to U, where Θ = Δ + UU'
function difference_precision_iht_grad_u!(grad_U, grad_Δ, A, Σx, Σy, U)
  diff_precision_grad_delta!(grad_Δ, A, Σx, Σy)
  rmul!(grad_Δ, 2)
  mul!(grad_U, grad_Δ, U)
end

# Δ ← Δ - η⋅G on the lower triangle
function _gradient_step!(Δ::SymmetricSparseIterate{T}, η, G::StridedMatrix{T}) where {T<:AbstractFloat}
  p = size(G, 1)
  for c=1:p, r=c:p
    Δ[r, c] = Δ[r, c] - T(η) * G[r, c]
  end
  Δ
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

# A = Σx⋅D⋅Σy without allocating, with tmp as storage
function _sandwich!(A::StridedMatrix{T}, tmp::StridedMatrix{T}, Σx::Symmetric{T}, D::StridedMatrix{T}, Σy::Symmetric{T}) where {T<:AbstractFloat}
  mul!(tmp, Σx, D)
  mul!(A, tmp, Σy)
end

# Θ = Δ + L as a dense matrix and A = Σx⋅Θ⋅Σy
function _latent_sandwich!(A, tmp, Θ, Σx, Δ, L, Σy)
  _densify!(Θ, Δ)
  Θ .+= L
  _sandwich!(A, tmp, Σx, Θ, Σy)
end

# the loss tr(A⋅Θ)/2 - tr(Θ⋅(Σy - Σx)) for a dense symmetric Θ with A = Σx⋅Θ⋅Σy
function _objective(A::StridedMatrix{T}, X::StridedMatrix{T}, Y::StridedMatrix{T}, Θ::StridedMatrix{T}) where {T<:AbstractFloat}
  v = zero(T)
  @inbounds for c in axes(Θ, 2), r in axes(Θ, 1)
    v += (A[r, c] / 2 + X[r, c] - Y[r, c]) * Θ[r, c]
  end
  v
end


"""
    difference_precision_iht(Σx, Σy, η, s; eps_tol=1e-5, maxiter=1000)

Direct estimation of the difference of precision matrices `Δ = Ωx - Ωy` with
at most `s` nonzero entries in its lower triangle, by projected gradient
descent on `tr(Σx⋅Δ⋅Σy⋅Δ)/2 - tr(Δ⋅(Σy - Σx))` with step size `η` and hard
thresholding. Stops when the objective changes by at most `eps_tol`.
Returns `(Δ, fvals)` with `Δ` a `SymmetricSparseIterate` and `fvals` the
objective values.
"""
function difference_precision_iht(
  Σx::Symmetric{T}, Σy::Symmetric{T},
  η::Real, s::Int;
  eps_tol=nothing, maxiter=nothing, epsTol=nothing, maxIter=nothing) where {T<:AbstractFloat}
  eps_tol = _naming_keyword(eps_tol, epsTol, 1e-5, :eps_tol, :epsTol)
  maxiter = _naming_keyword(maxiter, maxIter, 1000, :maxiter, :maxIter)

  size(Σx) == size(Σy) || throw(DimensionMismatch())

  p = size(Σx, 1)

  X, Y = Matrix(Σx), Matrix(Σy)
  gD = zeros(T, p, p)
  A = zeros(T, p, p)      # Σx⋅Δ⋅Σy
  D = zeros(T, p, p)      # Δ as a dense matrix
  tmp = zeros(T, p, p)
  Δ = SymmetricSparseIterate(T, p)

  fvals = T[]
  # gradient descent
  #
  fv = _objective(A, X, Y, D)
  push!(fvals, fv)
  for iter=1:maxiter
    # update Δ
    diff_precision_grad_delta!(gD, A, X, Y)
    _gradient_step!(Δ, η, gD)
    hard_threshold!(Δ, s)
    _sandwich!(A, tmp, Σx, _densify!(D, Δ), Σy)

    # check for convergence
    fv_new = _objective(A, X, Y, D)
    push!(fvals, fv_new)
    if abs(fv_new - fv) <= eps_tol
      break
    end
    fv = fv_new
  end

  (Δ, fvals)
end


####################################
#
# Direct difference estimation using iterative hard thresholding
#
#  now with low rank
#
####################################



"""
    difference_latent_precision_iht(Σx, Σy, ηΔ, ηU, s, r; options=IHTOptions(), callback=nothing)
    difference_latent_precision_iht!(Δ, L, U, gD, gU, Σx, Σy, ηΔ, ηU, s, r; options=IHTOptions(), callback=nothing)

Estimate the difference of precision matrices as `Δ + L`, with `Δ` sparse
(at most `s` nonzero entries in its lower triangle) and `L = UU'` of rank `r`,
by gradient descent on `tr(Σx⋅Θ⋅Σy⋅Θ)/2 - tr(Θ⋅(Σy - Σx))`, `Θ = Δ + UU'`,
with step sizes `ηΔ` and `ηU` and hard thresholding of `Δ`. The in-place
version starts from `Δ` and `U` and uses `gD` (`p × p`) and `gU` (`p × r`) as
gradient storage. `callback(Δ, L)`, if given, is called after every iteration.

Returns `(Δ, L, U, fvals, iter)`: the estimates, the objective values checked
every `options.check_every` iterations, and the number of iterations.
"""
function difference_latent_precision_iht!(
  Δ::SymmetricSparseIterate{T}, L::StridedMatrix{T},                 # these are pre-allocated
  U::StridedMatrix{T}, gD::StridedMatrix{T}, gU::StridedMatrix{T},   # these are pre-allocated
  Σx::Symmetric{T}, Σy::Symmetric{T},
  ηΔ, ηU, s::Int, r::Int;
  options::IHTOptions=IHTOptions(), callback=nothing) where {T<:AbstractFloat}

  eps_tol, maxiter, check_every = options.eps_tol, options.maxiter, options.check_every

  p = size(Σx, 1)
  X, Y = Matrix(Σx), Matrix(Σy)
  A = zeros(T, p, p)      # Σx⋅Θ⋅Σy
  Θ = zeros(T, p, p)      # Δ + L as a dense matrix
  tmp = zeros(T, p, p)
  mul!(L, U, U')
  _latent_sandwich!(A, tmp, Θ, Σx, Δ, L, Σy)

  fvals = T[]
  # gradient descent
  #
  fv = _objective(A, X, Y, Θ)
  push!(fvals, fv)
  iter = 0
  while iter < maxiter
    iter += 1
    # update Δ
    diff_precision_grad_delta!(gD, A, X, Y)
    max_gD = maximum(abs, gD)
    _gradient_step!(Δ, ηΔ, gD)
    hard_threshold!(Δ, s)

    # update A
    _latent_sandwich!(A, tmp, Θ, Σx, Δ, L, Σy)

    # update U and L
    difference_precision_iht_grad_u!(gU, gD, A, X, Y, U)
    @. U -= ηU * gU
    mul!(L, U, U')
    max_gU = maximum(abs, gU)

    # update A
    _latent_sandwich!(A, tmp, Θ, Σx, Δ, L, Σy)

    if callback !== nothing
      callback(Δ, L)
    end

    done = max(max_gD, max_gU) < eps_tol
    # check for convergence
    if mod(iter, check_every) == 0
      fv_new = _objective(A, X, Y, Θ)
      push!(fvals, fv_new)
      done = abs(fv_new - fv) <= eps_tol
      fv = fv_new
    end
    done && break
  end

  (Δ, L, U, fvals, iter)
end



# initialize Δ by hard thresholding inv(Σx) - inv(Σy) and U from the top r
# singular vectors of the remainder
function difference_latent_precision_iht_init(
  Σx::Symmetric{T}, Σy::Symmetric{T},
  s::Int, r::Int) where {T<:AbstractFloat}

  size(Σx) == size(Σy) || throw(DimensionMismatch())
  p = size(Σx, 1)

  gD = zeros(T, p, p)
  gU = zeros(T, p, r)
  Δ = SymmetricSparseIterate(T, p)

  # initialize Δ, L, U
  #
  tmp = Matrix(inv(Σx) - inv(Σy))
  hard_threshold!(Δ, tmp, s)
  for col_ind=1:p, row_ind=col_ind:p
    tmp[row_ind,col_ind] = tmp[row_ind,col_ind] - Δ[row_ind, col_ind]
    if row_ind != col_ind
      tmp[col_ind,row_ind] = tmp[row_ind,col_ind]
    end
  end
  F = svd(tmp)
  U = F.U[:,1:r] .* sqrt.(F.S[1:r])'
  L = U * U'

  (Δ, L, U, gD, gU)
end

function difference_latent_precision_iht(
  Σx::Symmetric{T}, Σy::Symmetric{T},
  ηΔ, ηU, s::Int, r::Int;
  options::IHTOptions=IHTOptions(), callback=nothing) where {T<:AbstractFloat}

  Δ, L, U, gD, gU = difference_latent_precision_iht_init(Σx, Σy, s, r)
  difference_latent_precision_iht!(Δ, L, U, gD, gU, Σx, Σy, ηΔ, ηU, s, r; options=options, callback=callback)
end






##
