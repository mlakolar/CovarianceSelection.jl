module RefitTest

using Test, Random, LinearAlgebra, SparseArrays
using CovarianceSelection
using ProximalBase: IndPattern, SymmetricSparseIterate

const rng = MersenneTwister(210)
const tight = ADMMOptions(; abstol=1e-12, reltol=1e-12, maxiter=10_000)
const methods = (ChordalCompletion(), DualCoordinate(; tol=1e-12),
  PrimalCoordinate(; tol=1e-12), RefitADMM(tight))

function pattern(p, edges)
  E = falses(p, p)
  for (i, j) in edges
    E[i, j] = E[j, i] = true
  end
  E
end

chain(p) = pattern(p, [(j, j + 1) for j in 1:p-1])
cycle(p) = pattern(p, [(j, mod1(j + 1, p)) for j in 1:p])

function grid(k)
  edges = Tuple{Int,Int}[]
  for j in 1:k, i in 1:k
    v = i + (j - 1) * k
    i < k && push!(edges, (v, v + 1))
    j < k && push!(edges, (v, v + k))
  end
  pattern(k^2, edges)
end

function fixture(kind, p)
  if kind == :chain
    cliques = [collect(j:j+1) for j in 1:p-1]
    separators = [[j] for j in 2:p-1]
  elseif kind == :band
    cliques = [collect(j:j+3) for j in 1:p-3]
    separators = [collect(j:j+2) for j in 2:p-3]
  elseif kind == :cliques
    cliques = [collect(j:min(j + 2, p)) for j in 1:3:p]
    separators = Vector{Int}[]
  else
    edges = [(kind == :star ? 1 : j ÷ 2, j) for j in 2:p]
    cliques = [[i, j] for (i, j) in edges]
    degrees = [count(c -> v in c, cliques) for v in 1:p]
    separators = [[v] for v in 1:p for _ in 1:degrees[v]-1]
  end
  E = falses(p, p)
  for c in cliques, j in c, i in c
    i == j || (E[i, j] = true)
  end
  if kind == :star
    # Duplicated leaves make the dual neighbourhood singular, while all edge
    # covariances remain positive definite and well conditioned.
    X = zeros(4, p)
    X[1, 1] = 1
    for j in 2:p
      X[1, j], X[2 + mod(j - 2, 3), j] = 0.2, 1
    end
    S = X'X
  else
    X = randn(rng, 4p, p)
    S = X'X / (4p)
  end
  S, E, cliques, separators
end

function oracle(S, cliques, separators)
  Ω = zeros(size(S))
  for (blocks, sign) in ((cliques, 1), (separators, -1)), c in blocks
    Ω[c, c] += sign * inv(Symmetric(S[c, c]))
  end
  Ω
end

function precision_error(Ω, expected, S)
  r = 1 ./ sqrt.(diag(S))
  maximum(abs.(Ω - expected) ./ (r * r'))
end

function certificate(Ω, S, E)
  r = 1 ./ sqrt.(diag(S))
  W = inv(Symmetric(Ω))
  maximum(abs(W[k, j] - S[k, j]) * r[k] * r[j]
    for j in axes(S, 2) for k in axes(S, 1) if k == j || E[k, j])
end

objective(Ω, S) = tr(S * Ω) - logdet(Symmetric(Ω))

function check_fit(S, E, method; expected=nothing, atol=1e-9, kkttol=1e-9, info=RefitInfo())
  estimate = covsel_refit(S, E, method; info)
  Ω = Matrix(estimate)
  @test estimate isa SymmetricSparseIterate{eltype(S)}
  @test issymmetric(Ω)
  @test isposdef(Symmetric(Ω))
  @test all(iszero(Ω[k, j]) for j in axes(E, 2) for k in axes(E, 1) if k != j && !E[k, j])
  @test info.converged
  expected === nothing || @test precision_error(Ω, expected, S) <= atol
  if method isa ChordalCompletion
    @test info.iterations == 0
    @test isnan(info.change) && isnan(info.kkt)
  else
    @test info.iterations >= 1
    @test isfinite(info.kkt)
    @test info.kkt <= kkttol
    @test info.kkt ≈ certificate(Ω, S, E) atol=100eps(eltype(S))
    if method isa RefitADMM
      @test isnan(info.change)
    else
      @test info.change <= max(eltype(S)(method.tol), 100eps(eltype(S)))
    end
  end
  Ω
end

@testset "refit" begin
  @testset "chordality" begin
    for kind in (:chain, :band, :cliques, :tree, :star)
      _, E, _, _ = fixture(kind, 10)
      before = copy(E)
      @test is_chordal(E)
      @test E == before
      perm = randperm(rng, 10)
      @test is_chordal(E[perm, perm])
    end
    @test is_chordal(trues(10, 10))
    @test is_chordal(falses(10, 10))
    @test is_chordal(falses(0, 0))
    @test !is_chordal(cycle(4))
    @test !is_chordal(grid(3))
  end

  @testset "clique and separator oracle" begin
    for p in (10, 30), kind in (:chain, :band, :cliques, :tree, :star)
      S, E, cliques, separators = fixture(kind, p)
      Ω = oracle(S, cliques, separators)
      perm = kind == :star ? [1; randperm(rng, p - 1) .+ 1] : randperm(rng, p)
      D = Diagonal(collect(range(0.8, 1.3; length=p)))
      for (S0, E0, Ω0) in ((S, E, Ω),
          (Matrix(Symmetric(D * S[perm, perm] * D)), E[perm, perm], inv(D) * Ω[perm, perm] * inv(D)))
        @test is_chordal(E0)
        for method in methods
          if kind == :star && method isa DualCoordinate
            @test_throws r"DualCoordinate broke down at node .*try PrimalCoordinate\(\)" covsel_refit(S0, E0, method)
          else
            check_fit(S0, E0, method; expected=Ω0)
          end
        end
      end
    end
  end

  @testset "non-chordal graphs and selected graphs" begin
    X = randn(rng, 60, 9)
    S = X'X / 60
    E = pattern(9, [(i, j) for j in 2:9 for i in 1:j-1 if rand(rng) < 0.25])
    E[1:4, 1:4] = cycle(4)
    fixtures = [(S, cycle(9)), (S, grid(3)), (S, E)]
    for n in (60, 7)
      X = randn(rng, n, 10)
      S0 = X'X / n
      Z = zeros(10, 10)
      covsel!(zeros(10, 10), Z, zeros(10, 10), S0, 0.25; options=tight, penalize_diag=false)
      push!(fixtures, (S0, IndPattern(Z).P))
    end
    population = inv(Symmetric(Matrix{Float64}(I, 10, 10) - 0.3 * cycle(10)))
    values, vectors = eigen(Symmetric(population))
    for n in (20, 9)
      X = if n > 10
        [sqrt(n) * cholesky(Symmetric(population)).U; zeros(n - 10, 10)]
      else
        sqrt(n) * Diagonal(sqrt.(values[2:end])) * vectors[:, 2:end]'
      end
      S0, Z = X'X / n, zeros(10, 10)
      covsel!(zeros(10, 10), Z, zeros(10, 10), S0, 0.1; options=tight, penalize_diag=false)
      E0 = IndPattern(Z).P
      @test !is_chordal(E0)
      push!(fixtures, (S0, E0))
    end
    admm = RefitADMM(ADMMOptions(; abstol=1e-10, reltol=1e-10, maxiter=10_000))
    for (index, (S0, E0)) in enumerate(fixtures)
      applicable = isposdef(Symmetric(S0)) ? (DualCoordinate(), PrimalCoordinate(), admm) : (PrimalCoordinate(), admm)
      fits = Matrix{Float64}[]
      for method in applicable
        info = RefitInfo()
        # ADMM stops on Frobenius residuals, which need not bound each scaled entry.
        kkttol = index == 5 && method isa RefitADMM ? 1e-8 : 1e-9
        push!(fits, check_fit(S0, E0, method; info, kkttol))
      end
      for Ω1 in fits, Ω2 in fits
        @test abs(objective(Ω1, S0) - objective(Ω2, S0)) <= 1e-8
      end
    end
  end

  @testset "singular three-node example" begin
    S = [1.0 0.5 0.5; 0.5 1 1; 0.5 1 1]
    E = pattern(3, [(1, 2), (1, 3)])
    expected = [5.0 -2 -2; -2 4 0; -2 0 4] / 3
    for method in (ChordalCompletion(), PrimalCoordinate(), RefitADMM(tight))
      check_fit(S, E, method; expected, atol=method isa RefitADMM ? 1e-8 : 1e-10)
    end
    @test_throws ErrorException("DualCoordinate broke down at node 1: the covariance iterate is not positive definite on its neighbourhood; the maximum likelihood estimate may still exist, try PrimalCoordinate()") covsel_refit(S, E, DualCoordinate())
  end

  @testset "empty, complete and one node" begin
    X = randn(rng, 20, 5)
    S = X'X / 20
    for method in methods
      check_fit(S, falses(5, 5), method; expected=Matrix(Diagonal(1 ./ diag(S))))
      check_fit(S, trues(5, 5), method; expected=inv(Symmetric(S)))
      for E in (falses(1, 1), trues(1, 1))
        check_fit(fill(2.0, 1, 1), E, method; expected=fill(0.5, 1, 1))
      end
    end
    @test RefitADMM().options isa ADMMOptions
    @test (DualCoordinate().tol, DualCoordinate().maxiter) == (1e-10, 1_000)
    @test (PrimalCoordinate().tol, PrimalCoordinate().maxiter) == (1e-10, 1_000)
    info = RefitInfo()
    @test (info.iterations, info.converged) == (0, false)
    @test isnan(info.change) && isnan(info.kkt)
  end

  @testset "graph forms and input preservation" begin
    S = [0.4^abs(i - j) for i in 1:6, j in 1:6]
    Z = zeros(6, 6)
    covsel!(zeros(6, 6), Z, zeros(6, 6), S, 0.2; options=tight)
    E = IndPattern(Z).P
    sparse_estimate = sparse(Z)
    sparse_estimate[1, 6] = sparse_estimate[6, 1] = 1
    sparse_estimate[1, 6] = sparse_estimate[6, 1] = 0
    @test nnz(sparse_estimate) > count(!iszero, sparse_estimate)
    @test !IndPattern(sparse_estimate).P[1, 6]
    for method in methods
      expected = Matrix(covsel_refit(S, E, method))
      for graph in (E, Matrix(E), Z, sparse_estimate, SymmetricSparseIterate(Z), IndPattern(E))
        before, sbefore = copy(graph isa IndPattern ? graph.P : graph), copy(S)
        Ω = Matrix(covsel_refit(S, graph, method))
        @test precision_error(Ω, expected, S) <= 1e-10
        @test S == sbefore
        @test (graph isa IndPattern ? graph.P : graph) == before
        @test is_chordal(graph) == is_chordal(E)
      end
      E0 = copy(E)
      E0[diagind(E0)] .= false
      @test Matrix(covsel_refit(S, E0, method)) ≈ expected atol=1e-10
      @test all(!E0[j, j] for j in 1:6)
      upper, lower = copy(S), copy(S)
      upper[6, 1], lower[1, 6] = NaN, Inf
      @test Matrix(covsel_refit(Symmetric(upper, :U), E, method)) ≈ expected atol=1e-10
      @test Matrix(covsel_refit(Symmetric(lower, :L), E, method)) ≈ expected atol=1e-10
      @test isnan(upper[6, 1]) && isinf(lower[1, 6])
    end
  end

  @testset "scales, permutations and Float32" begin
    S = [0.35^abs(i - j) for i in 1:6, j in 1:6]
    D = Diagonal([0.7, 1.3, 0.9, 1.5, 1.1, 0.8])
    perm = [4, 1, 6, 2, 5, 3]
    for E in (chain(6), cycle(6)), method in methods
      method isa ChordalCompletion && !is_chordal(E) && continue
      expected = Matrix(covsel_refit(S, E, method))
      scaled = Matrix(Symmetric(D * S * D))
      Ω = Matrix(covsel_refit(scaled, E, method))
      @test precision_error(Ω, inv(D) * expected * inv(D), scaled) <= 1e-10
      Ωp = Matrix(covsel_refit(S[perm, perm], E[perm, perm], method))
      @test precision_error(Ωp, expected[perm, perm], S[perm, perm]) <= 1e-10
      method32 = method isa RefitADMM ? RefitADMM(ADMMOptions(; abstol=1e-6, reltol=1e-6)) : method
      check_fit(Float32.(S), E, method32; expected, atol=1e-4, kkttol=1e-4)
    end
  end

  @testset "validation and failure messages" begin
    S, E = Matrix{Float64}(I, 3, 3), chain(3)
    for method in methods
      @test_throws DimensionMismatch("S must be a square matrix") covsel_refit(fill(NaN, 2, 3), E, method)
      @test_throws ArgumentError("S must not be empty") covsel_refit(zeros(0, 0), E, method)
      bad = copy(S)
      bad[1, 2] = NaN
      @test_throws ArgumentError("S must be finite") covsel_refit(bad, E, method)
      bad[1, 2] = Inf
      @test_throws ArgumentError("S must be finite") covsel_refit(bad, E, method)
      bad[1, 2], bad[1, 1] = 1.0, 0.0
      @test_throws ArgumentError("S must be symmetric") covsel_refit(bad, E, method)
      bad[1, 2] = 0.0
      @test_throws ArgumentError("the diagonal of S must be positive") covsel_refit(bad, falses(2, 2), method)
      bad[1, 1] = -1.0
      @test_throws ArgumentError("the diagonal of S must be positive") covsel_refit(bad, E, method)
      @test_throws DimensionMismatch("E must be a p × p pattern, the size of S") covsel_refit(S, falses(2, 2), method)
      @test_throws DimensionMismatch("E must be a p × p pattern, the size of S") covsel_refit(S, falses(3, 2), method)
      asym = copy(E)
      asym[1, 2] = false
      @test_throws ArgumentError("the pattern E must be symmetric") covsel_refit(S, asym, method)
    end
    @test_throws DimensionMismatch("E must be a square pattern") is_chordal(falses(2, 3))
    @test_throws DimensionMismatch("E must be a square pattern") is_chordal(IndPattern(falses(3)))
    @test_throws ArgumentError("the pattern E must be symmetric") is_chordal(Bool[0 1; 0 0])
    @test_throws MethodError covsel_refit(S, E)
    @test_throws ArgumentError("the graph E is not chordal; use DualCoordinate() or PrimalCoordinate()") covsel_refit(Matrix{Float64}(I, 4, 4), cycle(4), ChordalCompletion())
    @test_throws r"S is not positive definite on the clique \[1, 2\] of E; the maximum likelihood estimate does not exist" covsel_refit([1.0 1; 1 1], trues(2, 2), ChordalCompletion())
    @test_throws r"DualCoordinate broke down at node 1:.*try PrimalCoordinate\(\)" covsel_refit([1.0 1; 1 1], trues(2, 2), DualCoordinate())
    for constructor in (DualCoordinate, PrimalCoordinate)
      for tol in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError("tol must be positive and finite") constructor(; tol)
        @test_throws ArgumentError("tol must be positive and finite") constructor(tol, 1)
      end
      for maxiter in (0, -1)
        @test_throws ArgumentError("maxiter must be at least 1") constructor(; maxiter)
        @test_throws ArgumentError("maxiter must be at least 1") constructor(1e-10, maxiter)
      end
    end
  end

  @testset "non-convergence and reused info" begin
    S = [0.5^abs(i - j) for i in 1:9, j in 1:9]
    E = grid(3)
    for method in (PrimalCoordinate(; maxiter=1), DualCoordinate(; maxiter=1),
        RefitADMM(ADMMOptions(; maxiter=2)))
      info = RefitInfo(999, true, 0.0, 0.0)
      Ω = @test_logs (:warn, r"without converging") covsel_refit(S, E, method; info)
      @test !info.converged
      @test info.iterations == (method isa RefitADMM ? 2 : 1)
      @test isfinite(info.kkt) || info.kkt == Inf
      @test all(iszero(Ω[k, j]) for j in 1:9 for k in 1:9 if k != j && !E[k, j])
    end
    unequal = Matrix(Diagonal([8.92; fill(0.01, 8)]))
    info = RefitInfo()
    # Before convergence, over-relaxation can leave the projected iterate indefinite.
    early = RefitADMM(ADMMOptions(; maxiter=3, α=1.9, adaptive=false))
    @test_logs (:warn, r"without converging") covsel_refit(unequal, falses(9, 9), early; info)
    @test !info.converged && info.kkt == Inf && isnan(info.change)
    X = Float64[1 0 0 1 0 0; 0 1 0 0 1 0; 0 0 1 0 0 1]
    singular = X'X / 3
    @test_throws r"S is not positive definite on the clique \[.*\] of E; the maximum likelihood estimate does not exist" covsel_refit(singular, trues(6, 6), ChordalCompletion())
    info = RefitInfo()
    @test_logs (:warn, r"without converging") covsel_refit(singular, trues(6, 6), PrimalCoordinate(; maxiter=50); info)
    @test !info.converged && info.iterations == 50
    for method in methods
      covsel_refit(S, chain(9), method; info)
      @test info.converged
      if method isa ChordalCompletion
        @test info.iterations == 0 && isnan(info.change) && isnan(info.kkt)
      else
        @test info.iterations >= 1 && isfinite(info.kkt)
      end
    end
  end
end

end
