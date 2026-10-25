@testset "Naming compatibility" begin
  @test CovarianceSelection.fusedGraphicalLasso === fused_graphical_lasso
  @test CovarianceSelection.diffLoss === diff_loss
  @test CovarianceSelection.differencePrecisionIHT === difference_precision_iht
  canonical = IHTOptions(; eps_tol=1e-4, maxiter=7, check_every=2)
  legacy = IHTOptions(; epsTol=1e-4, maxIter=7, checkEvery=2)
  @test all(name -> hasproperty(legacy, name), (:epsTol, :maxIter, :checkEvery))
  @test fieldnames(IHTOptions) == (:eps_tol, :maxiter, :check_every)
  @test (canonical.eps_tol, canonical.maxiter, canonical.check_every) ==
        (legacy.epsTol, legacy.maxIter, legacy.checkEvery)
  @test_throws ArgumentError IHTOptions(; maxiter=7, maxIter=7)
  S = Symmetric([2.0 0.1; 0.1 1.0])
  T = Symmetric([1.5 0.2; 0.2 2.0])
  a = difference_precision_iht(S, T, 0.01, 1; eps_tol=1e-5, maxiter=3)
  b = CovarianceSelection.differencePrecisionIHT(S, T, 0.01, 1; epsTol=1e-5, maxIter=3)
  @test Matrix(a[1]) == Matrix(b[1])
  @test a[2] == b[2]
end
