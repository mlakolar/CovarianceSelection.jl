using Test,LinearAlgebra,Random
using CovarianceSelection
using CovarianceSelection.SingleIndexGraph
const MH_TEST=CovarianceSelection.SingleIndexGraph
@testset "MAVE selection with Stein coefficient refitting" begin
    rng=MersenneTwister(348151);Z=randn(rng,200,20)
    beta=zeros(20);beta[2:5]=[.5,.5,-.5,-.5]
    Y=(Z*beta)*[1. -2. .7];initial=normalize(beta+.02randn(rng,20))
    for refit_method in (:rrr,:ols)
        for method in (:threshold,:lasso,:adaptive_lasso)
            fit=select_mave_variables(Y,Z,initial;method,refit_method,
                   level=method==:threshold ? 1. : .1,iterative=method!=:threshold)
            @test fit.selected==[2,3,4,5]
            @test abs(dot(fit.direction,beta))>1-1e-10
            @test fit.status==:stein_refit && fit.iterations==0
            @test norm(fit.direction)≈1
        end
        # A new support orthogonal to the initial vector must not cause a
        # first-coordinate or orientation failure.
        old=zeros(20);old[1]=1.
        fit=MH_TEST._mave_stein_refit(Y,Z,[2,3,4,5],old,refit_method)
        @test abs(dot(fit.direction,beta))>1-1e-10
        @test fit.direction[1]==0.
    end
    @test_throws ArgumentError select_mave_variables(Y,Z,initial;refit_method=:unknown)
end
