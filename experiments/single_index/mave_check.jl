include("mave_common.jl")
using Test
include("../../test/mave.jl")
@testset "KernelSmoothing and paired Gaussian MAVE canary" begin
    for design in ("correlation","opposed_scale","fixed_diagonal","noncommuting")
        amplitude=design=="correlation" ? .55 : .6
        data=fg_data(810831,200,20,20,design,amplitude)
        Y=mp_responses(data.X,mp_bank(20,design=="opposed_scale" ? "mixed" : "complete",60,37151))
        b=mv_initial(Y,data.Z)
        t=data.Z*b;anchors=collect(1:200);h=.4
        localfit=MV._mave_smoothers(t,anchors,h)
        A,B=mv_smoother(t,anchors,h)
        @test A≈localfit.A atol=1e-10
        @test B≈localfit.B atol=1e-9
        @test maximum(abs.(vec(sum(A;dims=2)).-1))<1e-10
        @test maximum(abs.(A*t-t))<1e-10
        @test maximum(abs.(B*t.-1))<1e-9
        a=refine_probe_index(Y,data.Z,b;updates=2)
        c=refine_probe_index(Y,data.Z,b;updates=2,smoother=mv_smoother)
        @test abs(dot(a.direction,c.direction))>1-1e-10
        @test all(h->h.after<=h.before+1e-8,c.history)
    end
end
println("MAVE_CANARY_PASSED")
