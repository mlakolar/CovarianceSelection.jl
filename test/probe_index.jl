using Test, LinearAlgebra, Random
using CovarianceSelection.SingleIndexGraph

@testset "Multiple probe Gaussian index" begin
    rng=MersenneTwister(10491)
    Z=randn(rng,500,20)
    beta=zeros(20); beta[2:5].=0.5
    Y=(Z*beta)*[1.0 -2.0 0.5]
    for method in (:ols,:screen)
        fit=estimate_probe_index(Y,Z;method,support_size=6)
        @test abs(dot(fit.direction,beta))>1-1e-10
        @test norm(fit.direction)≈1
        @test fit.direction[1]≈0 atol=1e-12
    end
    u=estimate_probe_index(Y,Z).direction
    @test abs(dot(u,estimate_probe_index(Y .+ [4.0 7.0 -1.0],Z).direction))≈1
    @test abs(dot(u,estimate_probe_index(Y.*[2.0 0.1 8.0],Z).direction))≈1
    @test dot(u,estimate_probe_index(Y,Z;reference=-u).direction)≈-1
    @test_throws ArgumentError estimate_probe_index(Y,Z;method=:unknown)
    @test_throws DomainError estimate_probe_index(zeros(500,3),Z)
    @test_throws ArgumentError estimate_probe_index(Y,Z;method=:split,split=(1:250,251:500))
    @test_throws ArgumentError estimate_probe_index(Y,Z;method=:split,standardize=false,
                                                   split=(collect(1:250),collect(250:499)))
    for method in (:moment,:split,:second,:second_split)
        W=Y.+0.4*(Z*beta).^2*[1.0 0.7 -1.0]
        fit=estimate_probe_index(W,Z;method)
        @test all(isfinite,fit.direction)
        @test norm(fit.direction)≈1
    end
    H=probe_index_moments(Y,Z;order=2)
    @test all(H[:,:,l]≈H[:,:,l]' for l in axes(H,3))
    @test_throws ArgumentError probe_index_moments(Y,Z;order=3)
    @test_throws ArgumentError probe_index_moments(Y[1:499,:],Z)
    V=randn(rng,100000,4); b=[0.0,1.0,1.0,0.0]/sqrt(2)
    W=reshape((V*b).^2,:,1)
    HH=probe_index_moments(W,V;order=2)[:,:,1]
    @test opnorm(HH-2b*b')<0.12
    @test abs(dot(estimate_probe_index(W,V;method=:second).direction,b))>0.998
end
