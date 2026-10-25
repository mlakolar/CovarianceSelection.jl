using Test, LinearAlgebra, Random, CovarianceSelection
using CovarianceSelection.SingleIndexGraph
include("../experiments/single_index/gaussian_fixed_worker.jl")
@testset "Fixed-support Gaussian design and normal penalty" begin
    for p in (6,8,40),t in (-100.0,-0.5,0.0,0.5,100.0)
        P=fixed_precision_truth(t,p)
        @test isposdef(Symmetric(P))
        @test minimum(eigvals(Symmetric(P[1:p-2,1:p-2])))>=.27-1e-12
        @test (abs.(P).>0)==(abs.(fixed_precision_truth(0.0,p)).>0)
        @test .32-1e-12<=abs(P[1,2])<=.44+1e-12
        @test P[1,4]==0
        @test inv(Symmetric(P)) ≈ covariance_truth(t,p;graph="fixed_chain")
    end
    c=Dict("n"=>30000,"p"=>8,"d"=>6,"index_sparsity"=>4,"index_equal"=>true,
           "graph"=>"fixed_chain","design"=>"gaussian")
    D=generate_data(809010001,c)
    @test D.gamma[1:4]==fill(.5,4)
    @test D.t≈D.Z*D.gamma
    # Whiten each observation using an independently formed precision block.
    # This checks the vectorized generator, including its eigensystem orientation.
    whitened=zeros(30000,8)
    for i in 1:30000
        whitened[i,:]=cholesky(Symmetric(fixed_precision_truth(D.t[i],8))).U*D.X[i,:]
    end
    @test maximum(abs,whitened'*whitened/30000-I)<.035
    X=D.X[1:500,:];t=D.t[1:500]
    L=local_covariance(X,t,0.,1.2);A=[2,3,4,5,6,7,8]
    prob=.1;C=1.1
    lambda=GFSIG._regression_penalty(X,A,L,C,:normal,prob,1)
    # Express the same calibration directly in unnormalized kernel weights.
    w=biweight.(t./1.2)
    scale=maximum(sqrt(sum(w.^2 .* X[:,j].^2)/(sum(w)*sum(w.*X[:,j].^2))) for j in A)
    z=sqrt(2)*GFSIG.erfcinv(prob/length(A))
    @test lambda≈C*z*scale
    @test GFSIG._regression_penalty(X,A,L,2C,:normal,prob,1)≈2lambda
    f=fit_precision(X,t;points=[-.25,.25],entries=[(1,1),(1,2),(1,4)],
        bandwidth=.9,pilot_bandwidth=1.2,refit=true,penalty_rule=:normal,multiplier=C)
    @test all(isempty,f.failure)
    @test length(f.settings.nuisance_diagnostics)==6
    @test all(q->q.penalty>0 && q.effective_n>5,f.settings.nuisance_diagnostics)
    @test_throws ArgumentError fit_precision(X,t;points=[0.],entries=[(1,1)],bandwidth=.9,penalty_rule=:bad)
    @test_throws ArgumentError fit_precision(X,t;points=[0.],entries=[(1,1)],bandwidth=.9,penalty_rule=:normal,pilot=:aclime)
    @test_throws ArgumentError fit_precision(X,t;points=[0.],entries=[(1,1)],bandwidth=.9,penalty_probability=0.)
    # Exercise index and pointwise workers and independent count checks.
    root=mktempdir(); plan=Dict("stage"=>"index","seed_base"=>809000000,"shard_size"=>1,
        "penalty_probability"=>.1,"scenario"=>[merge(c,Dict("n"=>500,"id"=>1,"name"=>"check",
        "replicates"=>1,"index_methods"=>["mean_soft","mean_hard","mean_ols","mom_ols","oracle_support_ols"],
        "index_multipliers"=>[1.0]))])
    path=joinpath(root,"index.toml");open(io->TOML.print(io,plan),path,"w")
    gfshard(path,1,joinpath(root,"index"))
    @test length(readlines(joinpath(root,"index","shard-00001.csv")))==6
    plan["stage"]="pointwise";s=plan["scenario"][1]
    merge!(s,Dict("index_methods"=>["known"],"constants"=>[1.1],"h_factors"=>[1.0],"pilot_factors"=>[2.0],
          "variants"=>["plain","union_refit","corrected","oracle_nuisance"],"grid_min"=>-.25,"grid_max"=>.25,"grid_n"=>3,"draws"=>199))
    path=joinpath(root,"pointwise.toml");open(io->TOML.print(io,plan),path,"w")
    gfshard(path,1,joinpath(root,"pointwise"))
    pointrows=readlines(joinpath(root,"pointwise","shard-00001.csv"))
    @test length(pointrows)==65
    fields=[split(line,',') for line in pointrows[2:end]]
    @test all(r->length(r)==31,fields)
    @test all(r->r[18]=="\"\"",fields)
    unionrows=filter(r->r[8]=="\"union_refit\"" && r[9]=="\"pointwise\"",fields)
    @test length(unionrows)==9
    @test all(r->isfinite(parse(Float64,r[27])) && parse(Float64,r[27])>0,unionrows)
    # Test CSV consistency and aggregation with the actual production assessor.
    assessor=joinpath(@__DIR__,"../experiments/single_index/gaussian_fixed_assess.py")
    run(`python3 $assessor $(joinpath(root,"index.toml")) $(joinpath(root,"index"))`)
    run(`python3 $assessor $path $(joinpath(root,"pointwise"))`)
end
