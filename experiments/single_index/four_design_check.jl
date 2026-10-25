using Test
include("four_design_common.jl")
@testset "Four-design generation and aggregation" begin
    for family in ("correlation","opposed_scale","fixed_diagonal","noncommuting")
        a=family=="correlation" ? 0.55 : 0.6
        for d in (4,20)
            D=fg_data(453811+d,200,20,d,family,a)
            @test size(D.X)==(200,20)
            @test size(D.Z)==(200,d)
            @test norm(D.beta)==1
            @test count(!iszero,D.beta)==4
            d==20 && @test D.beta[1]==0
            Y=mp_responses(D.X,mp_bank(20,"mixed",60,480317))
            P=ag_prepare(Y,D.Z);W=Matrix{Float64}(I,size(Y,2),size(Y,2))
            for k in (d==4 ? (4,) : (10,20)),method in ("ols","rrr")
                u=ag_fit_matrix(P,W,method,k)
                @test abs(norm(u)-1)<1e-10
                @test count(!iszero,u)==k
            end
            @test abs(norm(ag_learned(Y,D.Z))-1)<1e-10
        end
    end
    # Known covariance factors used by the sampler must give inverse precision.
    G=fg_graphs(20)
    for family in ("fixed_diagonal","noncommuting"),t in (-1.,0.,1.)
        C=fg_curve(family,20,0.6,t,G);L=cholesky(Symmetric(C.O)).L
        F=inv(L');@test norm(F*F'-C.S)<1e-10
    end
end
