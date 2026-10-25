using CovarianceSelection
using Test
using Aqua

tests = [
  "utils",
  "covsel",
  "directDifference",
  "fusedDifference",
  "groupGraphicalLasso",
  "latent",
  "iht_solvers",
  "generate",
  "clime",
  "slassoinv",
  "single_index",
  "index_orientation",
  "probe_index",
  "mave",
  "mave_hybrid"
]

@testset "CovarianceSelection" begin
  for t in tests
    f = "$t.jl"
    println("* running $f ...")
    include(f)
  end

  @testset "Aqua" begin
    Aqua.test_all(CovarianceSelection)
  end
end
