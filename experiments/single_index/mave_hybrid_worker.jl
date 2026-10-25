# Feasible support selection followed by a first-order Stein coefficient refit.
# Frozen after development; confirmation observations/seeds match the primary run.
include("mave_common.jl")

function mh_refit(P,S;method=:rrr)
    d=size(P.Zc,2);isempty(S)&&error("empty selected support")
    A=P.Zc[:,S];root=cholesky(Symmetric(A'A/(size(A,1)-1))).L
    M=root\(A'P.V/(size(A,1)-1))
    B=method==:rrr ? M : root'\M
    v=ag_direction(B*B')
    method==:rrr && (v=root'\v)
    out=zeros(d);out[S]=normalize(v);out
end

function mh_worker(task)
    plan=TOML.parsefile("experiments/single_index/mave_confirmation.toml")
    tasks=[(s,r:min(r+plan["shard_size"]-1,s["replicates"])) for s in plan["scenario"]
           for r in 1:plan["shard_size"]:s["replicates"]]
    s,reps=tasks[task];out="hybrid";mkpath(out)
    path=joinpath(out,"part-$(lpad(task,4,'0')).csv");isfile(path)&&error("existing finalized shard")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,design,p,n,method,level,loss,missing,extra,all_true,exact_support,selection_status,selection_iterations,seconds,status,reason")
        for rep in reps
            seed=plan["seed"]+100000*s["id"]+rep;n=s["n"];p=s["p"]
            data=fg_data(seed,n,p,20,s["design"],s["amplitude"])
            Y=mp_responses(data.X,mp_bank(p,s["bank"],60,seed+47017))
            P=ag_prepare(Y,data.Z);b0=mv_initial(Y,data.Z)
            function emit(method,level,fitfun;selection_status="na",selection_iterations=0,shared_seconds=0.)
                started=time_ns();loss=1.;missing=-1;extra=-1;status="failed";reason=""
                try
                    b=fitfun();all(isfinite,b)&&abs(norm(b)-1)<1e-8||error("invalid fitted direction")
                    loss=1-clamp(abs(dot(b,data.beta)),0.,1.)^2;S=findall(!iszero,b)
                    missing=count(j->data.beta[j]!=0&&!(j in S),1:20)
                    extra=count(j->data.beta[j]==0&&j in S,1:20);status="ok"
                catch err
                    err isa InterruptException && rethrow()
                    reason=replace(sprint(showerror,err),','=>';','\n'=>' ')
                end
                println(io,join((s["id"],rep,s["design"],p,n,method,level,loss,missing,extra,Int(missing==0),Int(missing==0&&extra==0),selection_status,selection_iterations,shared_seconds+(time_ns()-started)/1e9,status,reason),','))
                flush(io)
            end
            emit("stein_screen10",0.,()->b0)
            for level in (.5,1.,1.5)
                S=findall(abs.(b0).>level*sqrt(log(20)/n))
                emit("stein_threshold_rrr_refit",level,()->mh_refit(P,S))
            end
            # One penalized quadratic update. Every variable is eligible.
            try
                started=time_ns();Q=MV._mave_system(P.V,data.Z,b0)
                weights=1 ./max.(abs.(b0),.05);lambda=.1maximum(abs.(Q.r)./weights)
                step=MV._mave_quadratic(Q.Q,Q.r;penalty=lambda.*weights)
                S=findall(!iszero,step.direction);elapsed=(time_ns()-started)/1e9
                emit("one_adaptive_select_rrr_refit",.1,()->mh_refit(P,S);
                      selection_status="one_update",selection_iterations=1,shared_seconds=elapsed)
            catch err
                emit("one_adaptive_select_rrr_refit",.1,()->throw(err))
            end
            # Identical penalized stage to the primary iterative adaptive MAVE arm.
            try
                started=time_ns()
                selected=refine_probe_index(Y,data.Z,b0;updates=60,bandwidth_factor=1.,
                                tolerance=1e-4,penalty_method=:adaptive_lasso,penalty_level=.1)
                elapsed=(time_ns()-started)/1e9
                for method in (:rrr,:ols)
                    emit("mave_select_$(method)_refit",.1,()->mh_refit(P,selected.selected;method);
                         selection_status=String(selected.status),selection_iterations=selected.iterations,shared_seconds=elapsed)
                end
            catch err
                for method in (:rrr,:ols);emit("mave_select_$(method)_refit",.1,()->throw(err));end
            end
            # Labeled diagnostic only. It never supplies support to a feasible arm.
            emit("oracle_support_rrr",0.,()->mh_refit(P,findall(!iszero,data.beta)))
        end
    end
    mv(path*".tmp",path)
end
if abspath(PROGRAM_FILE)==@__FILE__
    mh_worker(parse(Int,ARGS[1]))
end
