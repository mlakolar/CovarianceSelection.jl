include("four_design_common.jl")

function fg_worker(planpath,task,out)
    plan=TOML.parsefile(planpath)
    tasks=[(s,first:min(first+plan["shard_size"]-1,s["replicates"])) for s in plan["scenario"] for first in 1:plan["shard_size"]:s["replicates"]]
    s,reps=tasks[task];configs=filter(c->c["d"]==s["d"],plan["candidate"])
    mkpath(out);path=joinpath(out,"part-$(lpad(task,4,'0')).csv")
    isfile(path) && error("finalized shard exists")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,design,amplitude,p,d,n,bank,method,k,rho,loss,prediction,delta,bandwidth_ratio,weight_l1_small,weight_l1,weight_l1_large,missing,extra,seconds,status,reason")
        for rep in reps
            seed=plan["seed"]+100000*s["id"]+rep
            n,p,d=s["n"],s["p"],s["d"]
            data=fg_data(seed,n,p,d,s["design"],s["amplitude"])
            eval=randn(MersenneTwister(seed+71091),1000,d)
            function emit(bank,method,k,rho,fit)
                started=time_ns();status="ok";reason=""
                m=(loss=1.,prediction=2.,delta=sqrt(2.),ratio=sqrt(2.)/n^(-0.2),small=2.,l1=2.,large=2.,missing=-1,extra=-1)
                try
                    m=fg_metrics(fit(),data.beta,eval,n)
                catch err
                    err isa InterruptException && rethrow()
                    status="failed";reason=replace(sprint(showerror,err),','=>';','\n'=>' ')
                end
                secs=(time_ns()-started)/1e9
                println(io,join((s["id"],rep,s["design"],s["amplitude"],p,d,n,bank,method,k,rho,m.loss,m.prediction,m.delta,m.ratio,m.small,m.l1,m.large,m.missing,m.extra,secs,status,reason),','))
            end
            trace=vec(sum(abs2,data.X;dims=2))
            products=(vec(sum(data.X;dims=2)).^2 .-trace)./2
            oracle=vec(sum((data.X*data.score).*data.X;dims=2))
            for bank in ("global_trace","global_products","oracle_probe")
                y=bank=="global_trace" ? trace : bank=="global_products" ? products : oracle
                P=ag_prepare(reshape(y,:,1),data.Z)
                for c in filter(c->c["bank"]==bank,configs)
                    fit=c["method"]=="true_support_ols" ?
                        ()->ag_ols_scalar(y,data.Z;support=findall(!iszero,data.beta)) :
                        ()->ag_fit_matrix(P,ones(1,1),"ols",c["k"])
                    emit(bank,c["method"],c["k"],c["rho"],fit)
                end
            end
            for bank in ("projections","mixed","complete")
                Y=mp_responses(data.X,mp_bank(p,bank,60,seed+47017))
                started=time_ns();P=ag_prepare(Y,data.Z);weights=Dict(ag_weights(P));setup=(time_ns()-started)/1e9
                L=size(Y,2);identity=Matrix{Float64}(I,L,L)
                for c in filter(c->c["bank"]==bank,configs)
                    method=c["method"];rho=c["rho"];k=c["k"]
                    fit=if startswith(method,"weighted_")
                        ()->ag_fit_matrix(P,weights[rho],replace(method,"weighted_"=>""),k)
                    elseif method=="learned_crossfit"
                        ()->ag_learned(Y,data.Z;rho)
                    elseif method in ("split","second")
                        ()->estimate_probe_index(Y,data.Z;method=Symbol(method)).direction
                    else
                        ()->ag_fit_matrix(P,identity,method,k)
                    end
                    emit(bank,method,k,rho,fit)
                end
                println(io,join((s["id"],rep,s["design"],s["amplitude"],p,d,n,bank,"setup",0,0.,0.,0.,0.,0.,0.,0.,0.,0,0,setup,"timing",""),','))
            end
        end
    end
    mv(path*".tmp",path)
end

if abspath(PROGRAM_FILE)==@__FILE__
    fg_worker(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
