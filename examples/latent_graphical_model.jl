# A chain graph confounded by one latent factor connected to all observed
# variables. Marginalizing it subtracts a positive semidefinite rank one
# matrix from the sparse precision: Ω = S - L. The graphical lasso of the
# observed variables then sees spurious edges from the latent factor.
#
# Run from the repository root (see examples/README.md for the setup):
#
#   julia --project=examples examples/latent_graphical_model.jl

using CovarianceSelection, ProximalSolvers
using LinearAlgebra, Printf, Random, Statistics

Random.seed!(1)

p, n = 20, 2000
S = Matrix(SymTridiagonal(ones(p), fill(0.4, p - 1)))
u = ones(p) / sqrt(p)
L = u * u'
Ω = S - L
X = CovarianceSelection.generate_data(Ω, n)
Ŝ = cov(X)
@printf("p = %d, n = %d, chain edges = %d, L[i,j] = %.2f, λmin(Ω) = %.3f\n\n",
  p, n, p - 1, L[1, 1], eigmin(Symmetric(Ω)))

# Fit the sparse conditional graph and the latent component together, with
# backtracking steps from the default diagonal S and small deterministic U.
Ŝ_fit, L̂, U, fvals, iter = covsel_latent_iht!(Ŝ, p - 1, 1;
  method = MonotoneAccProxGradDescent(),
  options = ProximalOptions(; maxiter = 5000, ftol = 1e-12, xtol = 1e-12, gradtol = 1e-10))

# Count each edge once, in the lower off-diagonal triangle. The true edges
# join adjacent nodes of the chain.
function edge_counts(Θ)
  p = size(Θ, 1)
  tp, fp = 0, 0
  for c=1:p-1, r=c+1:p
    if abs(Θ[r, c]) > 1e-4
      if r == c + 1
        tp += 1
      else
        fp += 1
      end
    end
  end
  tp, fp
end

tp, fp = edge_counts(Ŝ_fit)
û = U[:, 1] / norm(U)
@printf("latent graphical model (%d iterations): %d true positive edges, %d false positive edges\n",
  iter, tp, fp)
@printf("  |cos ∠(û, u)| = %.3f, ‖Ŝ_fit - S‖_F = %.3f, ‖L̂ - L‖_F = %.3f\n\n",
  abs(dot(û, u)), norm(Matrix(Ŝ_fit) - S), norm(L̂ - L))

# Compare with the graphical lasso at the closest edge count to p - 1.
# Penalties descend, so the first minimum breaks ties toward the larger λ.
# The diagonal is penalized, the default of covselpath.
λs = exp10.(range(-0.5, -2.5; length = 30))
path = covselpath(Ŝ, λs)
counts = edge_counts.(path)
best = argmin([abs(tp + fp - (p - 1)) for (tp, fp) in counts])
tp, fp = counts[best]
@printf("graphical lasso (λ = %.4f): %d true positive edges, %d false positive edges\n",
  λs[best], tp, fp)
@printf("  selected edge count = %d, target = %d\n", tp + fp, p - 1)
if best == 1 && minimum(tp + fp for (tp, fp) in counts) > p - 1
  @printf("  The graphical lasso does not reach %d edges on this grid.\n", p - 1)
end
