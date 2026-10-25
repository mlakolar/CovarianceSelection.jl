include("mave_common.jl")

function mv_worker(planpath,task,out)
    plan=TOML.parsefile(planpath)
    tasks=[(s,first:min(first+plan["shard_size"]-1,s["replicates"]))
           for s in plan["scenario"] for first in 1:plan["shard_size"]:s["replicates"]]
    s,reps=tasks[task];mkpath(out);path=joinpath(out,"part-$(lpad(task,4,'0')).csv")
    isfile(path) && error("finalized shard exists")
    trajectory=joinpath(out,"trajectory-$(lpad(task,4,'0')).csv")
    tio=open(trajectory*".tmp","w")
    println(tio,"scenario,replicate,design,p,n,bandwidth_factor,iteration,loss,change")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,design,amplitude,p,d,n,bank,method,bandwidth_factor,level,loss,prediction,delta,missing,extra,exact_support,all_true,iterations,status,seconds,change,validation,selection_status,selection_iterations,reason")
        for rep in reps
            seed=plan["seed"]+100000*s["id"]+rep
            n,p,d=s["n"],s["p"],20
            data=fg_data(seed,n,p,d,s["design"],s["amplitude"])
            Y=mp_responses(data.X,mp_bank(p,s["bank"],60,seed+47017))
            # Truth is used only below in the reporting function, never in a fit.
            function emit(method,h,level,fitfun;validation=NaN)
                started=time_ns();reason="";b=nothing;fit=nothing
                status="failed";iterations=0;change=NaN
                loss=1.;prediction=2.;delta=sqrt(2.);missing=-1;extra=-1
                try
                    fit=fitfun();b=fit.direction
                    c=clamp(abs(dot(b,data.beta)),0.,1.)
                    loss=1-c*c;prediction=2-2c;delta=sqrt(2-2c)
                    selected=findall(!iszero,b)
                    missing=count(j->data.beta[j]!=0&&!(j in selected),1:d)
                    extra=count(j->data.beta[j]==0&&j in selected,1:d)
                    status=String(fit.status);iterations=fit.iterations
                    hasproperty(fit,:history) && !isempty(fit.history) && (change=fit.history[end].change)
                    if method=="full_mave"
                        for step in fit.history
                            e=1-clamp(abs(dot(step.direction,data.beta)),0.,1.)^2
                            println(tio,join((s["id"],rep,s["design"],p,n,h,step.iteration,e,step.change),','))
                        end
                        flush(tio)
                    end
                catch err
                    err isa InterruptException && rethrow()
                    fit=nothing;status="failed";loss=1.;prediction=2.;delta=sqrt(2.);missing=-1;extra=-1
                    reason=replace(sprint(showerror,err),','=>';','\n'=>' ')
                end
                seconds=(time_ns()-started)/1e9
                println(io,join((s["id"],rep,s["design"],s["amplitude"],p,d,n,s["bank"],method,h,level,
                                loss,prediction,delta,missing,extra,Int(missing==0&&extra==0),Int(missing==0),
                                iterations,status,seconds,change,validation,
                                fit!==nothing&&hasproperty(fit,:selection_status) ? String(fit.selection_status) : "na",
                                fit!==nothing&&hasproperty(fit,:selection_iterations) ? fit.selection_iterations : 0,reason),','))
                flush(io)
                fit
            end
            init=emit("stein_screen10",0.,0.,()->(;direction=mv_initial(Y,data.Z),status=:stein,iterations=0))
            emit("stein_dense20",0.,0.,()->(;direction=mv_initial(Y,data.Z,d),status=:stein,iterations=0))
            for h in plan["bandwidth_factors"]
                emit("one_update",h,0.,()->refine_probe_index(Y,data.Z,init.direction;updates=1,bandwidth_factor=h))
                dense=emit("full_mave",h,0.,()->refine_probe_index(Y,data.Z,init.direction;
                                    updates=get(plan,"dense_updates",plan["max_updates"]),bandwidth_factor=h,tolerance=plan["tolerance"]))
                # Fixed penalties: not selected by truth or cross-validation.
                for method in ("threshold","lasso","adaptive_lasso")
                    levels=method=="threshold" ? plan["threshold_levels"] : plan["lasso_levels"]
                    for level in levels
                        emit(method*"_refit",h,level,()->select_mave_variables(Y,data.Z,dense.direction;
                             method=Symbol(method),level,refit_updates=plan["max_updates"],bandwidth_factor=h,
                             tolerance=plan["tolerance"]))
                    end
                end
                if get(plan,"iterative_selection",false) && h==plan["holdout_bandwidth"]
                    for rule in ("lasso","adaptive_lasso")
                        emit("iterative_"*rule*"_refit",h,.1,()->select_mave_variables(Y,data.Z,init.direction;
                             method=Symbol(rule),level=.1,iterative=true,refit_updates=plan["max_updates"],
                             bandwidth_factor=h,tolerance=plan["tolerance"]))
                    end
                end
                if plan["holdout"] && h==plan["holdout_bandwidth"]
                    started=time_ns()
                    try
                        selection=mv_holdout(Y,data.Z;seed=seed+91731,bandwidth_factor=h,updates=plan["max_updates"])
                        elapsed=(time_ns()-started)/1e9
                        for name in ("holdout_min","holdout_1se")
                            chosen=selection.results[name]
                            # Shared selection work is recorded once in a separate timing row.
                            emit(name,h,Float64(chosen.k),()->(;chosen.fit...,selection_status=Symbol(chosen.training_status),selection_iterations=0);validation=chosen.validation)
                        end
                        println(io,join((s["id"],rep,s["design"],s["amplitude"],p,d,n,s["bank"],"holdout_setup",h,0.,
                                         0.,0.,0.,-1,-1,0,0,0,"timing",elapsed,NaN,NaN,"na",0,""),','))
                    catch err
                        err isa InterruptException && rethrow()
                        for name in ("holdout_min","holdout_1se")
                            emit(name,h,0.,()->throw(err))
                        end
                    end
                end
            end
        end
    end
    close(tio)
    mv(trajectory*".tmp",trajectory)
    mv(path*".tmp",path)
end

if abspath(PROGRAM_FILE)==@__FILE__
    mv_worker(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
