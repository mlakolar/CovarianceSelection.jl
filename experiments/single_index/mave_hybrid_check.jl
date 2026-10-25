include("mave_hybrid_worker.jl")
using Test
@testset "Stein refits on data-selected supports" begin
    rng=MersenneTwister(349151);Z=randn(rng,300,20)
    beta=zeros(20);beta[2:5]=[.5,.5,-.5,-.5]
    Y=(Z*beta)*[1. -2. .7];P=ag_prepare(Y,Z);S=findall(!iszero,beta)
    for method in (:rrr,:ols)
        b=mh_refit(P,S;method)
        @test abs(dot(b,beta))>1-1e-10
        @test all(iszero,b[setdiff(1:20,S)])
        @test norm(b)≈1
    end
    Y .+=randn(rng,size(Y)...);P=ag_prepare(Y,Z)
    S=sort(sortperm(vec(sum(abs2,P.B;dims=2));rev=true)[1:6])
    for method in (:rrr,:ols)
        b=mh_refit(P,S;method)
        reference=ag_fit_matrix(P,Matrix{Float64}(I,3,3),String(method),6)
        @test abs(dot(b,reference))>1-1e-10
    end
    @test_throws ErrorException mh_refit(P,Int[])
end
