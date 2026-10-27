@testset "Naming compatibility" begin
  @test CovarianceSelection.fusedGraphicalLasso === fused_graphical_lasso
  @test CovarianceSelection.diffLoss === diff_loss
  @test CovarianceSelection.differencePrecisionIHT === difference_precision_iht
  @test CovarianceSelection.differenceLatentPrecisionIHT === difference_latent_precision_iht
  @test CovarianceSelection.differenceLatentPrecisionIHT! === difference_latent_precision_iht!
  @test CovarianceSelection.covselLatentIHT! === covsel_latent_iht!
  @test !isdefined(CovarianceSelection, :IHTOptions)
  S = Symmetric([2.0 0.1; 0.1 1.0])
  T = Symmetric([1.5 0.2; 0.2 2.0])
  @test_throws MethodError difference_precision_iht(S, T, 0.01, 1)
  @test_throws MethodError difference_latent_precision_iht(S, T, 0.01, 0.01, 1, 1)
  @test_throws MethodError covsel_latent_iht!(S, 0.05, 0.05, 1, 1)
end
