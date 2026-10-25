using Test,CovarianceSelection,LinearAlgebra
include("../experiments/single_index/dgp.jl")
@testset "Fixed-support hub and known-index controls" begin
    for t in (-100.,-.5,0.,.5,100.)
        P=fixed_hub_precision(t,16;degree=6)
        @test isposdef(Symmetric(P))
        @test minimum(eigvals(Symmetric(P[1:14,1:14])))>=1.4-.44sqrt(6)-1e-12
        @test P[2,3]==0
        @test count(!iszero,P[1,:])-1==6
        @test (P.!=0)==(fixed_hub_precision(0.,16;degree=6).!=0)
    end
    c=Dict("n"=>20000,"p"=>16,"d"=>20,"index_sparsity"=>4,"index_equal"=>true,"graph"=>"fixed_hub","graph_degree"=>6)
    D=generate_data(907101001,c);W=zeros(20000,16)
    for i in 1:20000
        W[i,:]=cholesky(Symmetric(fixed_hub_precision(D.t[i],16;degree=6))).U*D.X[i,:]
    end
    @test maximum(abs,W'W/20000-I)<.04
end
