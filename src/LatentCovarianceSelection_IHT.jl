"""
    IHTOptions(; eps_tol=1e-2, maxiter=500, check_every=10)

Options for the iterative hard thresholding solvers. The iterations stop when
the largest gradient entry falls below `eps_tol`, or when the objective, checked
every `check_every` iterations, changes by at most `eps_tol`.
"""
struct IHTOptions
 eps_tol::Float64
 maxiter::Int64
 check_every::Int64
end


function _naming_keyword(value, legacy, default, canonical::Symbol, old::Symbol)
  if legacy !== nothing
    value === nothing || throw(ArgumentError("supply only $canonical, not both $canonical and $old"))
    Base.depwarn("keyword $old is deprecated; use $canonical", old)
    return legacy
  end
  value === nothing ? default : value
end

function IHTOptions(; eps_tol::Union{Nothing,Real}=nothing,
                    maxiter::Union{Nothing,Integer}=nothing,
                    check_every::Union{Nothing,Integer}=nothing,
                    epsTol::Union{Nothing,Real}=nothing,
                    maxIter::Union{Nothing,Integer}=nothing,
                    checkEvery::Union{Nothing,Integer}=nothing)
  IHTOptions(_naming_keyword(eps_tol, epsTol, 1e-2, :eps_tol, :epsTol),
             _naming_keyword(maxiter, maxIter, 500, :maxiter, :maxIter),
             _naming_keyword(check_every, checkEvery, 10, :check_every, :checkEvery))
end

function Base.getproperty(options::IHTOptions, name::Symbol)
  canonical = name === :epsTol ? :eps_tol : name === :maxIter ? :maxiter :
              name === :checkEvery ? :check_every : name
  canonical === name || Base.depwarn("field $name is deprecated; use $canonical", name)
  getfield(options, canonical)
end
Base.propertynames(::IHTOptions, private::Bool=false) =
  (fieldnames(IHTOptions)..., :epsTol, :maxIter, :checkEvery)



  # gradient of tr(Σ(S + L)) - log det(S + L) with respect to S
  function covsel_latent_iht_grad_s!(gS, Σ, S, L)
    gS .= Symmetric(Σ) - inv(Symmetric(S + L))
    gS
  end

  # gradient with respect to U, where L = UU'
  function covsel_latent_iht_grad_u!(gU, gS, Σ, S, L, U)
    covsel_latent_iht_grad_s!(gS, Σ, S, L)
    rmul!(gS, 2.)
    mul!(gU, gS, U)
  end

  # minimize tr(Σ(S + L)) - log det(S + L) over S with at most s nonzero
  # off-diagonal entries (lower triangle) and L = UU' of rank r, by gradient
  # descent with step sizes ηS and ηU and hard thresholding of S
  function covsel_latent_iht!(Σ::StridedMatrix, ηS, ηU, s, r;
    eps_tol=nothing, maxiter=nothing, check_every=nothing, callback=nothing,
    epsTol=nothing, maxIter=nothing, checkEvery=nothing)
    eps_tol = _naming_keyword(eps_tol, epsTol, 1e-2, :eps_tol, :epsTol)
    maxiter = _naming_keyword(maxiter, maxIter, 100, :maxiter, :maxIter)
    check_every = _naming_keyword(check_every, checkEvery, 5, :check_every, :checkEvery)

    p = size(Σ, 1)

    gS = zeros(p, p)
    gU = zeros(p, r)

    S = zeros(p, p)

    # init
    tmp = Matrix(inv(Symmetric(Σ)))
    hard_threshold!(S, tmp, s, false)

    tmp .= tmp - S
    F = svd(tmp)
    U = F.U[:,1:r] .* sqrt.(F.S[1:r])'
    L = U * U'


    fvals = Float64[]
    # gradient descent
    #
    fv = covsel_objective(Σ, S, L)
    push!(fvals, fv)
    iter = 0
    while iter < maxiter
      iter += 1
      # update S
      covsel_latent_iht_grad_s!(gS, Σ, S, L)
      @. S -= ηS * gS
      hard_threshold!(S, S, s, false)
      max_gS = maximum( abs.(gS) )

      # update U and L
      covsel_latent_iht_grad_u!(gU, gS, Σ, S, L, U)
      @. U -= ηU * gU
      mul!(L, U, U')
      max_gU = maximum( abs.(gU) )
      if callback !== nothing
        callback(Σ, S, L)
      end

      done = max(max_gS, max_gU) < eps_tol
      # check for convergence
      if mod(iter, check_every) == 0
        fv_new = covsel_objective(Σ, S, L)
        push!(fvals, fv_new)
        done = abs(fv_new - fv) <= eps_tol
        fv = fv_new
      end
      done && break
    end

    (S, L, U, fvals, iter)
  end
