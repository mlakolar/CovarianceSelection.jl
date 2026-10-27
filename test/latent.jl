module LatentTest

using Test, Random, LinearAlgebra
import CovarianceSelection, ProximalBase
using ProximalBase: SymmetricSparseIterate, nnz
using ProximalSolvers: ProximalOptions, MonotoneAccProxGradDescent

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

    opt = ProximalOptions(; maxiter = 5000, ftol = 1e-14, xtol = 1e-14, gradtol = 1e-12)
    calls = Ref(0)
    Δ, L, U, fvals, iter = CovarianceSelection.difference_latent_precision_iht(Symmetric(Σx), Symmetric(Σy), 1, r;
        options = opt, callback = (Δ, L) -> (calls[] += 1))
    @test calls[] == iter
    @test 1 <= iter <= opt.maxiter
    @test L ≈ U * U'
    @test size(U) == (p, r)
    @test nnz(Δ) <= 1

    # the fit improves on the starting point and the residual of the
    # estimating equation is small
    Δ0, U0 = SymmetricSparseIterate(p), zeros(p, r)
    U0[1, 1] = 0.1
    res0 = CovarianceSelection.diff_loss(Symmetric(Σx), Δ0, U0, Symmetric(Σy), 2)
    res = CovarianceSelection.diff_loss(Symmetric(Σx), Δ, U, Symmetric(Σy), 2)
    @test res < res0
    @test fvals[end] <= fvals[1]
    @test Matrix(Δ) + L ≈ Θ atol=0.05

    # the in-place version continues from the given iterates
    Δ1, L1, U1, _, iter1 = CovarianceSelection.difference_latent_precision_iht!(copy(Δ), copy(L), copy(U),
        Symmetric(Σx), Symmetric(Σy), 1;
        options = ProximalOptions(; maxiter = 5, ftol = 1e-14, xtol = 1e-14, gradtol = 1e-12))
    # automatic stopping may finish before the iteration limit
    @test 1 <= iter1 <= 5
    @test CovarianceSelection.diff_loss(Symmetric(Σx), Δ1, U1, Symmetric(Σy), 2) <= res + 1e-8

    # the objective values are those of the loss at Θ = Δ + UU'
    Θ = Matrix(Δ1) + U1 * U1'
    @test CovarianceSelection.difference_latent_precision_iht!(copy(Δ1), zeros(p, p), copy(U1),
        Symmetric(Σx), Symmetric(Σy), 1; options = ProximalOptions(; maxiter = 2, callback = _ -> true))[4][1] ≈
      tr(Σx * Θ * Σy * Θ) / 2 - tr(Θ * (Σy - Σx))

    # only the triangle that each Symmetric refers to is read
    Mx, My = copy(Σx), copy(Σy)
    Mx[2, 1] = My[1, 3] = 1e3
    short = ProximalOptions(; maxiter = 20, ftol = 1e-14, xtol = 1e-14, gradtol = 1e-12)
    a = CovarianceSelection.difference_latent_precision_iht(Symmetric(Σx), Symmetric(Σy), 1, r; options = short)
    b = CovarianceSelection.difference_latent_precision_iht(Symmetric(Mx, :U), Symmetric(My, :L), 1, r; options = short)
    @test Matrix(b[1]) ≈ Matrix(a[1])
    @test b[2] ≈ a[2]
  end

  @testset "covsel_latent_iht!" begin
    p = 10
    S = Matrix(SymTridiagonal(ones(p), fill(0.4, p - 1)))
    u = ones(p) / sqrt(p)
    L = u * u'
    Σ = inv(Symmetric(S - L))

    calls = Ref(0)
    positive_definite = Bool[]
    callback = (Σ, S, L) -> begin
      calls[] += 1
      push!(positive_definite, isposdef(Symmetric(Matrix(S) - L)))
    end
    Ŝ, L̂, U, fvals, iter = CovarianceSelection.covsel_latent_iht!(Σ, p - 1, 1;
        method = MonotoneAccProxGradDescent(),
        options = ProximalOptions(; maxiter = 5000, ftol = 1e-12, xtol = 1e-12, gradtol = 1e-10),
        callback)
    @test calls[] == iter
    @test 1 <= iter <= 5000
    @test all(positive_definite)
    @test all(diff(fvals) .<= 1e-10 .* (1 .+ abs.(fvals[1:end-1])))
    @test L̂ ≈ U * U'
    @test CovarianceSelection.covsel_objective(Σ, Ŝ, L̂) ≈ fvals[end]
    @test issymmetric(Ŝ)
    @test isposdef(Symmetric(Matrix(Ŝ) - L̂))
    support = Set((i, j) for j in 1:p for i in j+1:p if !iszero(Ŝ[i, j]))
    @test support == Set((i + 1, i) for i in 1:p-1)
    @test abs(dot(U[:, 1] / norm(U), u)) >= 0.95
    @test CovarianceSelection.covsel_objective(Σ, Ŝ, L̂) >= CovarianceSelection.covsel_objective(Σ, inv(Σ)) - 1e-8
  end
end

end
