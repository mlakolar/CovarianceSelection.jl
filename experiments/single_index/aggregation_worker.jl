include("aggregation_common.jl")

function ag_diagnostics(u,beta,Zeval,n)
    all(isfinite,u) && abs(norm(u)-1)<1e-7 || error("invalid direction")
    dot(u,beta)<0 && (u=-u)
    delta=norm(u-beta); t=Zeval*beta;that=Zeval*u
    errs=Float64[]
    for c in (0.75,1.0,1.5)
        h=c*n^(-0.2)
        w=exp.(-0.5.*(t./h).^2);w ./=sum(w)
        v=exp.(-0.5.*(that./h).^2);v ./=sum(v)
        push!(errs,sum(abs,w-v))
    end
    m=mp_metrics(u,beta)
    (;loss=m.loss,prediction=m.prediction,delta,ratio=delta/n^(-0.2),l1small=errs[1],l1=errs[2],l1large=errs[3])
end

function ag_worker(planpath,task,out)
    plan=TOML.parsefile(planpath)
    tasks=[(s,first:min(first+plan["shard_size"]-1,s["replicates"]))
      for s in plan["scenario"] for first in 1:plan["shard_size"]:s["replicates"]]
    s,reps=tasks[task];mkpath(out)
    path=joinpath(out,"part-$(lpad(task,4,'0')).csv")
    isfile(path) && error("finalized shard exists")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,n,p,design,signal,beta_kind,bank,method,k,rho,loss,prediction,delta,bandwidth_ratio,weight_l1_small,weight_l1,weight_l1_large,seconds,status,reason")
        for rep in reps
            seed=plan["seed"]+100000*s["id"]+rep
            n,p=s["n"],s["p"]
            D=mp_data(MersenneTwister(seed),n,p,s["design"];signal=s["signal"],beta_kind=s["beta_kind"])
            Zeval=randn(MersenneTwister(seed+91961),1000,20)
            function emit(bank,method,k,rho,fit)
                started=time_ns();status="ok";reason=""
                m=(loss=1.0,prediction=2.0,delta=sqrt(2.),ratio=NaN,l1small=2.,l1=2.,l1large=2.)
                try
                    u=fit();m=ag_diagnostics(u,D.beta,Zeval,n)
                catch e
                    e isa InterruptException && rethrow()
                    status="failed";reason=replace(sprint(showerror,e),','=>';','\n'=>' ')
                end
                secs=(time_ns()-started)/1e9
                println(io,join((s["id"],rep,n,p,s["design"],s["signal"],s["beta_kind"],bank,method,k,rho,
                    m.loss,m.prediction,m.delta,m.ratio,m.l1small,m.l1,m.l1large,secs,status,reason),','))
            end
            q=ag_oracle_probe(D.X,s["design"],s["signal"])
            emit("oracle_probe","ols",20,0.,()->ag_ols_scalar(q,D.Z))
            emit("oracle_probe","true_support_ols",4,0.,()->ag_ols_scalar(q,D.Z;support=findall(!iszero,D.beta)))
            banks=p<=40 ? ("projections","mixed","complete") : ("projections","mixed")
            for bank in banks
                Y=mp_responses(D.X,mp_bank(p,bank,60,seed+47017))
                started=time_ns();P=ag_prepare(Y,D.Z);L=size(Y,2)
                setupsecs=(time_ns()-started)/1e9
                identity=Matrix{Float64}(I,L,L)
                for method in ("ols","moment","rrr"),k in (10,20)
                    emit(bank,method,k,1.,()->ag_fit_matrix(P,identity,method,k))
                end
                # A shared eigendecomposition makes the three fixed levels cheap.
                # Setup timing is reported separately to avoid pretending it is free.
                started=time_ns();weights=ag_weights(P);weightsecs=(time_ns()-started)/1e9
                for (rho,W) in weights,method in ("ols","rrr"),k in (10,20)
                    emit(bank,"weighted_"*method,k,rho,()->ag_fit_matrix(P,W,method,k))
                end
                emit(bank,"learned_crossfit",20,0.2,()->ag_learned(Y,D.Z;rho=0.2))
                println(io,join((s["id"],rep,n,p,s["design"],s["signal"],s["beta_kind"],bank,"setup",0,0.,
                    0.,0.,0.,0.,0.,0.,0.,setupsecs+weightsecs,"timing",""),','))
            end
        end
    end
    mv(path*".tmp",path)
end

if abspath(PROGRAM_FILE)==@__FILE__
    ag_worker(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
