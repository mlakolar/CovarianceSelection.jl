"""
    ChordalCompletion()

Maximum likelihood on a chordal graph by a product of regressions in a perfect
elimination ordering. Each regression uses a clique `Cᵥ` formed by a vertex and
its later neighbours. The cost is Σᵥ |Cᵥ|³, without iterations or a final
inverse. It throws an `ArgumentError` if the graph is not chordal. It also
throws if `S` is not positive definite on a clique, since the maximum
likelihood estimate then does not exist. There are no options or stopping
tolerances.
"""
struct ChordalCompletion <: RefitMethod end

"""
    is_chordal(E)

Test whether the undirected graph `E` is chordal by maximum cardinality search
and a perfect elimination ordering test, in O(p²) time on the Boolean pattern.
Accepts the graph forms of [`covsel_refit`](@ref), including `IndPattern`;
the pattern must be square and symmetric. Its diagonal is ignored, and the
input is not modified. The empty graph is chordal.
"""
is_chordal(E) = _refit_order(_refit_pattern(E)) !== nothing

function _refit_order(P)
  p = size(P, 1)
  weights, order, position = zeros(Int, p), zeros(Int, p), zeros(Int, p)
  for t in p:-1:1
    v, best = 0, -1
    for k in 1:p
      if position[k] == 0 && weights[k] > best
        v, best = k, weights[k]
      end
    end
    order[t], position[v] = v, t
    for k in 1:p
      P[k, v] && position[k] == 0 && (weights[k] += 1)
    end
  end
  for (t, v) in enumerate(order)
    follower, first = 0, p + 1
    for k in 1:p
      if P[k, v] && t < position[k] < first
        follower, first = k, position[k]
      end
    end
    follower == 0 && continue
    for k in 1:p
      if P[k, v] && position[k] > t && k != follower && !P[k, follower]
        return nothing
      end
    end
  end
  order
end

function _refit_clique_error(c)
  throw(ArgumentError("S is not positive definite on the clique $(sort(c)) of E; the maximum likelihood estimate does not exist"))
end

function _refit(S::Matrix{T}, P, ::ChordalCompletion, info) where {T}
  order = _refit_order(P)
  order === nothing && throw(ArgumentError("the graph E is not chordal; use DualCoordinate() or PrimalCoordinate()"))
  p = size(S, 1)
  position = invperm(order)
  Ω = zeros(T, p, p)
  for (t, v) in enumerate(order)
    N = [k for k in 1:p if P[k, v] && position[k] > t]
    β = S[N, v]
    if !isempty(N)
      factor = cholesky(Symmetric(S[N, N]); check=false)
      isposdef(factor) || _refit_clique_error([v; N])
      ldiv!(factor, β)
    end
    σ² = S[v, v] - dot(view(S, v, N), β)
    σ² > zero(T) || _refit_clique_error([v; N])
    Ω[v, v] += one(T) / σ²
    for (a, k) in enumerate(N)
      Ω[k, v] -= β[a] / σ²
      Ω[v, k] = Ω[k, v]
      for b in a:length(N)
        l = N[b]
        Ω[l, k] += β[a] * β[b] / σ²
        Ω[k, l] = Ω[l, k]
      end
    end
  end
  info.iterations, info.converged, info.change, info.kkt = 0, true, NaN, NaN
  Ω
end
