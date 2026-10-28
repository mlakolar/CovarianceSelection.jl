using CovarianceSelection
using Test
using Aqua
using LinearAlgebra

tests = [
  "utils",
  "covsel",
  "refit",
  "directDifference",
  "fusedDifference",
  "multigroupFusedGraphicalLasso",
  "danaher2014",
  "group_graphical_lasso",
  "latent",
  "iht_solvers",
  "generate",
  "clime",
  "slassoinv"
]

@testset "CovarianceSelection" begin
  for t in tests
    f = "$t.jl"
    println("* running $f ...")
    include(f)
  end

  include("naming.jl")

  @testset "Aqua" begin
    Aqua.test_all(CovarianceSelection)
  end
end
