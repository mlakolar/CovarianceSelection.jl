# Bounded algorithm diagnostic, not a retuned primary comparison.
include("mave_common.jl")
using Printf
plan=TOML.parsefile("experiments/single_index/mave_development.toml")
cases=[(s,r) for s in plan["scenario"] if s["p"]==40&&s["n"] in (400,800) for r in 1:2]
function convergence_task(task)
    s,rep=cases[task];n=s["n"];seed=plan["seed"]+100000*s["id"]+rep
    data=fg_data(seed,n,40,20,s["design"],s["amplitude"])
    Y=mp_responses(data.X,mp_bank(40,s["bank"],60,seed+47017))
    initial=mv_initial(Y,data.Z);started=time_ns()
    fit=refine_probe_index(Y,data.Z,initial;updates=2000,bandwidth_factor=1.,tolerance=1e-4)
    loss(b)=1-clamp(abs(dot(b,data.beta)),0.,1.)^2
    h=fit.history;N=length(h)
    path=joinpath("convergence","part-$(lpad(task,3,'0')).csv")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,design,n,iterations,status,stein_loss,loss1,loss20,loss60,loss200,loss500,loss1000,loss2000,change,cycle2,cycle5,seconds")
        errors=[loss(h[min(i,N)].direction) for i in (1,20,60,200,500,1000,2000)]
        cycle2=N>=3 ? norm(h[end].direction-h[end-2].direction) : NaN
        cycle5=N>=6 ? norm(h[end].direction-h[end-5].direction) : NaN
        println(io,join((s["id"],rep,s["design"],n,N,fit.status,loss(initial),errors...,
                         h[end].change,cycle2,cycle5,(time_ns()-started)/1e9),','))
    end
    mv(path*".tmp",path)
end
mkpath("convergence")
convergence_task(parse(Int,ARGS[1]))
