using Test, LinearAlgebra, Random
using CovarianceSelection.SingleIndexGraph
const MAVE_TEST=CovarianceSelection.SingleIndexGraph

@testset "Projected multivariate MAVE" begin
    rng=MersenneTwister(810387);Z=randn(rng,100,6)
    beta=normalize([0.,1.,-1.,1.,0.,0.]);Y=(Z*beta)*[1. -2. .7]
    initial=normalize([0.1,1.,-.6,0.,.2,.1])
    V=(Y .-sum(Y;dims=1)/100)./sqrt.(sum(abs2,Y.-sum(Y;dims=1)/100;dims=1)/99)
    P=MAVE_TEST._mave_system(V,Z,initial)
    Q=zeros(6,6);r=zeros(6)
    # Independent direct accumulation of the stacked multivariate regression.
    for (ii,i) in enumerate(P.anchors),j in 1:100
        dz=Z[j,:]-Z[i,:];s=P.slope[ii,:];a=P.intercept[ii,:]
        w=P.W[ii,j]/(length(P.anchors)*size(Y,2))
        Q .+=w*dot(s,s).*(dz*dz')
        r .+=w*dot(s,V[j,:]-a).*dz
    end
    @test Q≈P.Q atol=1e-10
    @test r≈P.r atol=1e-10
    step=MAVE_TEST._mave_quadratic(P.Q,P.r)
    @test norm(P.Q*step.direction-P.r)<1e-8
    penalty=fill(0.05maximum(abs,P.r),6)
    sparse=MAVE_TEST._mave_quadratic(P.Q,P.r;penalty)
    gradient=P.Q*sparse.direction-P.r
    @test maximum(iszero(sparse.direction[j]) ? max(abs(gradient[j])-penalty[j],0.) :
                  abs(gradient[j]+penalty[j]*sign(sparse.direction[j])) for j in 1:6)<1e-7
    fit=refine_probe_index(Y,Z,initial;updates=100)
    @test abs(dot(beta,fit.direction))>0.999
    @test abs(fit.direction[4])>0.3 # Missing initial coordinate can reenter.
    @test abs(fit.direction[1])<0.01 # First coordinate is truly zero.
    @test norm(fit.direction)≈1
    @test fit.status==:converged
    @test all(x->x.after<=x.before+1e-8,fit.history)
    neg=refine_probe_index(Y,Z,-initial;updates=100)
    @test abs(dot(fit.direction,neg.direction))>1-1e-10
    one=refine_probe_index(Y,Z,initial;updates=1)
    @test one.iterations==1 && one.status==:one_update
    restricted=refine_probe_index(Y,Z,initial;updates=10,support=[2,3,4])
    @test all(iszero,restricted.direction[[1,5,6]])
    for method in (:lasso,:adaptive_lasso,:threshold)
        selected=select_mave_variables(Y,Z,fit.direction;method,
                                      level=method==:threshold ? 1. : 0.1)
        @test selected.selected==[2,3,4]
        @test abs(dot(beta,selected.direction))>0.999
    end
    for method in (:lasso,:adaptive_lasso)
        iter=select_mave_variables(Y,Z,fit.direction;method,level=.1,iterative=true)
        @test iter.selected==[2,3,4]
        @test abs(dot(beta,iter.direction))>0.999
        @test iter.selection_status==:converged
    end
    @test_throws ArgumentError refine_probe_index(Y,Z,zeros(6))
    @test_throws ArgumentError refine_probe_index(Y,Z,initial;bandwidth_factor=0.)
    @test_throws ArgumentError refine_probe_index(Y,Z,initial;support=Int[])
    @test_throws ArgumentError select_mave_variables(Y,Z,fit.direction;method=:unknown)
end
