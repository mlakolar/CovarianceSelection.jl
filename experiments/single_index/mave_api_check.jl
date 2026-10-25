include("mave_hybrid_worker.jl")
using Test
@testset "Public hybrid agrees with frozen experiment" begin
    data=fg_data(3537512,400,40,20,"correlation",.55)
    Y=mp_responses(data.X,mp_bank(40,"complete",60,3584529))
    P=ag_prepare(Y,data.Z);initial=mv_initial(Y,data.Z)
    pen=refine_probe_index(Y,data.Z,initial;updates=60,penalty_method=:adaptive_lasso,penalty_level=.1)
    for method in (:rrr,:ols)
        fit=select_mave_variables(Y,data.Z,initial;method=:adaptive_lasso,level=.1,
                                 iterative=true,refit_updates=60,refit_method=method)
        @test fit.selected==pen.selected
        @test fit.selection_status==pen.status
        @test abs(dot(fit.direction,mh_refit(P,pen.selected;method)))>1-1e-10
    end
end
