module LatentTest

using Test, Random, LinearAlgebra
import CovarianceSelection, ProximalBase
using ProximalBase: SymmetricSparseIterate, nnz

Random.seed!(11)

@testset "latent" begin

  @testset "difference_latent_precision_iht" begin
    p, r = 8, 1
    Σx = [0.5^abs(i-j) for i=1:p, j=1:p]
    S = zeros(p, p)
    S[1, 2] = S[2, 1] = 0.4
    u = randn(p); u ./= norm(u)
    Θ = S + 0.3 * u * u'
    Σy = Matrix(inv(Symmetric(inv(Σx) - Θ)))
    @test isposdef(Symmetric(Σy))

    η = 0.5 / (opnorm(Σx) * opnorm(Σy))
    opt = CovarianceSelection.IHTOptions(;eps_tol=1e-10, maxiter=20000, check_every=10)
    calls = Ref(0)
    Δ, L, U, fvals, iter = CovarianceSelection.difference_latent_precision_iht(Symmetric(Σx), Symmetric(Σy), η, η, 1, r;
        options=opt, callback=(Δ, L) -> (calls[] += 1))
    @test calls[] == iter
    @test 1 <= iter <= opt.maxiter
    @test L ≈ U * U'
    @test size(U) == (p, r)
    @test nnz(Δ) <= 1

    # the fit improves on the starting point and the residual of the
    # estimating equation is small
    Δ0, L0, U0, _, _ = CovarianceSelection.difference_latent_precision_iht_init(Symmetric(Σx), Symmetric(Σy), 1, r)
    res0 = CovarianceSelection.diff_loss(Symmetric(Σx), Δ0, U0, Symmetric(Σy), 2)
    res = CovarianceSelection.diff_loss(Symmetric(Σx), Δ, U, Symmetric(Σy), 2)
    @test res < res0
    @test fvals[end] <= fvals[1]
    @test Matrix(Δ) + L ≈ Θ atol=0.05

    # the in-place version continues from the given iterates
    Δ1, L1, U1, _, iter1 = CovarianceSelection.difference_latent_precision_iht!(copy(Δ), copy(L), copy(U), zeros(p, p), zeros(p, r),
        Symmetric(Σx), Symmetric(Σy), η, η, 1, r; options=CovarianceSelection.IHTOptions(;eps_tol=1e-10, maxiter=5, check_every=10))
    @test iter1 == 5
    @test CovarianceSelection.diff_loss(Symmetric(Σx), Δ1, U1, Symmetric(Σy), 2) <= res + 1e-8

    # the objective values are those of the loss at Θ = Δ + UU'
    Θ = Matrix(Δ1) + U1 * U1'
    @test CovarianceSelection.difference_latent_precision_iht!(copy(Δ1), zeros(p, p), copy(U1), zeros(p, p), zeros(p, r),
        Symmetric(Σx), Symmetric(Σy), η, η, 1, r; options=CovarianceSelection.IHTOptions(;maxiter=1))[4][1] ≈
      tr(Σx * Θ * Σy * Θ) / 2 - tr(Θ * (Σy - Σx))

    # only the triangle that each Symmetric refers to is read
    Mx, My = copy(Σx), copy(Σy)
    Mx[2, 1] = My[1, 3] = 1e3
    short = CovarianceSelection.IHTOptions(;eps_tol=1e-10, maxiter=20, check_every=10)
    a = CovarianceSelection.difference_latent_precision_iht(Symmetric(Σx), Symmetric(Σy), η, η, 1, r; options=short)
    b = CovarianceSelection.difference_latent_precision_iht(Symmetric(Mx, :U), Symmetric(My, :L), η, η, 1, r; options=short)
    @test Matrix(b[1]) ≈ Matrix(a[1])
    @test b[2] ≈ a[2]
  end

  @testset "covsel_latent_iht!" begin
    p, r = 6, 1
    Ω = Matrix(2.0I, p, p)
    Ω[1, 2] = Ω[2, 1] = 0.5
    u = ones(p) / sqrt(p)
    Σ = inv(Ω + 0.5 * u * u')

    calls = Ref(0)
    S, L, U, fvals, iter = CovarianceSelection.covsel_latent_iht!(Σ, 0.05, 0.05, 1, r;
        eps_tol=1e-10, maxiter=200, check_every=5, callback=(Σ, S, L) -> (calls[] += 1))
    @test calls[] == iter
    @test L ≈ U * U'
    @test issymmetric(S)
    @test isposdef(Symmetric(S + L))
    @test count(!iszero, tril(S, -1)) <= 1
    @test fvals[end] <= fvals[1]
    @test CovarianceSelection.covsel_objective(Σ, S, L) >= CovarianceSelection.covsel_objective(Σ, inv(Σ)) - 1e-8
  end
end

end
