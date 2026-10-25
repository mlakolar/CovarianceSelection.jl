include("mave_common.jl")
using Test
@testset "Finite validation candidate selection" begin
    best,sparse,eligible=mv_support_choice([NaN,1.1,1.],[NaN,.1,.1])
    @test best==3 && sparse==2 && eligible==[false,true,true]
    best,sparse,eligible=mv_support_choice([1.05,1.,Inf],[.1,.1,Inf])
    @test best==2 && sparse==1 && eligible==[true,true,false]
    @test_throws DomainError mv_support_choice([NaN,Inf],[NaN,Inf])
end
# Reproduce the one failed validation draw, with the same primary inputs.
seed=2937511+100000*8+179
D=fg_data(seed,1600,40,20,"correlation",.55)
Y=mp_responses(D.X,mp_bank(40,"complete",60,seed+47017))
fit=mv_holdout(Y,D.Z;seed=seed+91731,bandwidth_factor=1.,updates=60)
mkpath("holdout_guard")
open("holdout_guard/candidates.csv","w") do io
    println(io,"size,score,standard_error,nonfinite_losses,eligible")
    for j in eachindex(fit.sizes)
        println(io,join((fit.sizes[j],fit.scores[j],fit.standard_errors[j],fit.nonfinite_losses[j],fit.eligible[j]),','))
    end
end
open("holdout_guard/repaired-case.csv","w") do io
    println(io,"scenario,replicate,method,selected_size,loss,status,scope")
    for name in ("holdout_min","holdout_1se")
        r=fit.results[name];loss=1-clamp(abs(dot(r.fit.direction,D.beta)),0.,1.)^2
        println(io,join((8,179,name,r.k,loss,r.fit.status,"diagnostic_only_primary_failures_retained"),','))
    end
end
println("FINITE_VALIDATION_GUARD_CHECKED")
