include("multiprobe_common.jl")

function mp_worker(planpath,task,outdir)
    plan=TOML.parsefile(planpath)
    tasks=[(s,first:min(first+plan["shard_size"]-1,s["replicates"]))
           for s in plan["scenario"] for first in 1:plan["shard_size"]:s["replicates"]]
    scenario,reps=tasks[task]
    configs=if haskey(plan,"candidate")
        [(id=c["id"],bank=c["bank"],L=c["L"],method=c["method"],k=c["k"]) for c in plan["candidate"]]
    else
        mp_configs()
    end
    mkpath(outdir)
    path=joinpath(outdir,"part-$(lpad(task,4,'0')).csv")
    isfile(path) && error("refusing to overwrite completed shard: $path")
    tmp=path*".tmp"
    open(tmp,"w") do io
        println(io,"scenario,replicate,n,p,design,signal,beta_kind,even,method,bank,L,k,selector,selected_id,loss,prediction,angle,missing,extra,seconds,status,reason")
        for rep in reps
            sid=scenario["id"]
            seed=plan["seed"]+100000sid+rep
            n,p=scenario["n"],scenario["p"]
            design=scenario["design"]; signal=get(scenario,"signal",1.0)
            bk=get(scenario,"beta_kind","equal"); even=get(scenario,"even",false)
            data=mp_data(MersenneTwister(seed),n,p,design;signal,beta_kind=bk,even)
            probe_seed=seed+47017
            cache=Dict{Tuple{String,Int},Matrix{Float64}}()
            losses=fill(Inf,length(configs)); records=Vector{Any}(undef,length(configs))
            for (i,c) in enumerate(configs)
                started=time_ns()
                reason=""; status="ok"; fit=nothing
                metrics=(loss=NaN,prediction=NaN,angle=NaN,missing=-1,extra=-1)
                try
                    Y=get!(cache,(c.bank,c.L)) do
                        mp_responses(data.X,mp_bank(p,c.bank,c.L,probe_seed))
                    end
                    # Compare fitting time on cached probes for every method.
                    started=time_ns()
                    fit=mp_fit(Y,data.Z,c)
                    metrics=mp_metrics(fit.direction,data.beta)
                    losses[i]=metrics.loss
                catch err
                    err isa InterruptException && rethrow()
                    status="failed";reason=replace(sprint(showerror,err),','=>';','\n'=>' ')
                    metrics=(loss=NaN,prediction=NaN,angle=NaN,missing=-1,extra=-1)
                end
                seconds=(time_ns()-started)/1e9
                records[i]=(c,metrics,seconds,status,reason)
                mp_write(io,scenario,rep,n,p,design,signal,bk,even,c,"fixed",c.id,metrics,seconds,status,reason)
            end
            if get(plan,"validate",false)
                started=time_ns(); selected=0
                try
                    # Validation and full-data fits use the same probe seed.
                    selected,scores=mp_validate(data.X,data.Z,configs,probe_seed+11)
                    c,m,secs,status,reason=records[selected]
                    mp_write(io,scenario,rep,n,p,design,signal,bk,even,c,"validation",c.id,m,
                             (time_ns()-started)/1e9,status,reason)
                catch err
                    err isa InterruptException && rethrow()
                    c=configs[1]; m=(loss=NaN,prediction=NaN,angle=NaN,missing=-1,extra=-1)
                    reason=replace(sprint(showerror,err),','=>';','\n'=>' ')
                    mp_write(io,scenario,rep,n,p,design,signal,bk,even,c,"validation","none",m,
                             (time_ns()-started)/1e9,"failed",reason)
                end
                c,m,secs,status,reason=records[argmin(losses)]
                mp_write(io,scenario,rep,n,p,design,signal,bk,even,c,"oracle",c.id,m,0.0,status,reason)
            end
        end
    end
    mv(tmp,path)
end

function mp_write(io,s,r,n,p,design,signal,bk,even,c,selector,selected,m,seconds,status,reason)
    println(io,join((s["id"],r,n,p,design,signal,bk,even,c.method,c.bank,c.L,c.k,selector,
                    selected,m.loss,m.prediction,m.angle,m.missing,m.extra,seconds,status,reason),','))
end

if abspath(PROGRAM_FILE)==@__FILE__
    mp_worker(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
