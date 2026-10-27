# A difference of two precision matrices that is sparse plus low rank, as when
# one group has a latent factor that the other lacks. It is estimated by
# iterative hard thresholding (difference_latent_precision_iht) and compared with
# the sparse direct difference estimator.
#
# Run from the repository root (see examples/README.md for the setup):
#
#   julia --project=examples examples/latent_difference.jl

using CovarianceSelection, ProximalBase, ProximalSolvers
using LinearAlgebra, Printf, Random, Statistics
using UnicodePlots

Random.seed!(3)

# Ωx is a chain graph. In group y two edges change and a latent variable is
# marginalized out, which subtracts a positive semidefinite rank one matrix L
# from the precision matrix. So Δ = Ωx - Ωy = S + L with S sparse.
p = 30
Ωx = Matrix(1.0I, p, p)
for i=1:p-1
  Ωx[i, i+1] = Ωx[i+1, i] = 0.3
end
S = zeros(p, p)
S[4, 3] = S[3, 4] = 0.3         # edge (3, 4) is absent in group y
S[20, 10] = S[10, 20] = -0.3    # edge (10, 20) is present only in group y
L = Matrix(CovarianceSelection.gen_latent_model(Ωx, 1, 0.2)[2])   # random, rank one, ‖L‖ = 0.2 λmax(Ωx)
Ωy = Ωx - S - L
@assert isposdef(Symmetric(Ωy))

n = 5000
Sx = Symmetric(cov(CovarianceSelection.generate_data(Ωx, n)))
Sy = Symmetric(cov(CovarianceSelection.generate_data(Ωy, n)))
@printf("p = %d, n = %d per group, ‖S‖_F = %.3f, ‖L‖_F = %.3f, largest |L[i,j]| = %.3f\n\n",
  p, n, norm(S), norm(L), maximum(abs, L))

# difference_latent_precision_iht minimizes the direct difference loss
#
#   ℓ(Θ) = tr(Sx Θ Sy Θ) / 2 - tr(Θ (Sy - Sx)),   Θ = Δ + UU',
#
# by gradient steps on Δ and U with a common step size chosen by backtracking.
# Each step projects Δ onto at most s lower-triangle coordinates, including
# the diagonal. The monotone accelerated method keeps the loss nonincreasing.
# The sparsity s and the rank r are tuning parameters; here they are the truth.
s, r = 2, 1

# The callback observes accepted iterates and records ‖Θₖ - Θₖ₋₁‖_F.
# Initially Δ = 0 and U has 0.1 on its first r diagonal entries.
steps = Float64[]
Θprev = Ref(zeros(p, p))
for i=1:r
  Θprev[][i, i] = 0.1^2
end
function recordStep(Δ, L)
  Θ = Matrix(Δ) + L
  push!(steps, norm(Θ - Θprev[]))
  Θprev[] = Θ
end
Δ, L̂, U, fvals, iter = difference_latent_precision_iht(Sx, Sy, s, r;
  method = MonotoneAccProxGradDescent(),
  options = ProximalOptions(; maxiter = 10_000, ftol = 1e-12, xtol = 1e-12, gradtol = 1e-10),
  callback = recordStep)

support = sort([(r, c) for c=1:p for r=c:p if Δ[r, c] != 0])
u, û = eigvecs(Symmetric(L))[:, end], U[:, 1] / norm(U)
println("sparse plus low rank, ", iter, " iterations")
println("  support of Δ̂: ", support)
@printf("  ‖Δ̂ - S‖_F = %.3f, ‖L̂ - L‖_F = %.3f, |cos ∠(û, u)| = %.3f\n", norm(Matrix(Δ) - S), norm(L̂ - L), abs(dot(û, u)))
@printf("  ‖Δ̂ + L̂ - (S + L)‖_F = %.3f\n\n", norm(Matrix(Δ) + L̂ - S - L))

println(lineplot(1:length(steps), log10.(steps); xlabel="iteration", ylabel="log10 step",
  title="Convergence", width=50, height=10))

# The sparse estimator has no low rank part, so it either leaves L out or
# spreads it over many small entries.
λ = 2 * sqrt(log(p) / n)
Δsparse = difference_precision_active_shooting(Sx, Sy, NormL1(λ))
@printf("\nsparse only (λ = %.3f): %d nonzero entries in the lower triangle, ‖Δ̂ - (S + L)‖_F = %.3f\n",
  λ, nnz(Δsparse), norm(Matrix(Δsparse) - S - L))
