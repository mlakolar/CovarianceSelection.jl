using Test
include("aggregation_common.jl")
rng=MersenneTwister(324981)
@testset "Aggregation identities" begin
    n=500;d=20;L=12
    Z=randn(rng,n,d);beta=mp_beta("mixed");a=randn(rng,L)
    Y=(Z*beta)*a'
    P=ag_prepare(Y,Z)
    for method in ("ols","moment","rrr")
        u=ag_fit_matrix(P,Matrix{Float64}(I,L,L),method,d)
        method!="moment" && @test abs(dot(u,beta))>1-1e-10
        @test abs(norm(u)-1)<1e-12
    end
    for (rho,W) in ag_weights(P)
        @test minimum(eigvals(Symmetric(W)))>0
        u=ag_fit_matrix(P,W,"rrr",20)
        @test abs(dot(u,beta))>1-1e-10
    end
    # RRR must equal the generalized eigenproblem, not plain OLS/SVD.
    Yn=Y+randn(rng,n,L);P=ag_prepare(Yn,Z)
    v=eigen(Symmetric(P.M*P.M'),Symmetric(P.Zc'P.Zc/(n-1))).vectors[:,end]
    u=ag_fit_matrix(P,Matrix{Float64}(I,L,L),"rrr",d)
    @test abs(dot(normalize(v),u))>1-1e-10
    Q=Matrix(qr(randn(rng,L,L)).Q)
    P2=merge(P,(N=P.N*Q,B=P.B*Q,M=P.M*Q,V=P.V*Q))
    @test abs(dot(u,ag_fit_matrix(P2,Matrix{Float64}(I,L,L),"rrr",d)))>1-1e-10
    us=ag_fit_matrix(P,Matrix{Float64}(I,L,L),"rrr",10)
    @test count(!iszero,us)==10
    @test abs(norm(ag_learned(Yn,Z))-1)<1e-10
    @test abs(dot(ag_learned(Y,Z),beta))>0.9
end
