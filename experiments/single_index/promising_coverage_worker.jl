# Paired pointwise coverage with independent index training. No truth enters fitting.
include("four_design_common.jl")

function pc_targets(design,p)
    if design=="opposed_scale"
        return [(1,1),(3,3),(1,2),(1,3)],["diagonal_plus","diagonal_minus","edge","zero"]
    elseif design=="correlation"
        return [(1,1),(2,2),(1,2),(1,3)],["diagonal_1","diagonal_2","edge","zero"]
    end
    g=fg_graphs(p)
    edge=findfirst(>(0),g.B)
    varying=minmax(Tuple(edge)...)
    zero=first((j,k) for j in 1:p for k in j+1:p if g.A[j,k]==0 && g.B[j,k]==0)
    [(1,1),varying,(1,2),zero],["diagonal","varying_edge","constant_edge","zero"]
end

function pc_index(train,perm,seed,bank,k)
    original=randperm(MersenneTwister(seed+39017),size(train.X,2))
    aligned=train.X[:,invperm(original)[perm]]
    Y=mp_responses(aligned,mp_bank(size(aligned,2),bank,60,seed+47017))
    P=ag_prepare(Y,train.Z)
    b=ag_fit_matrix(P,Matrix{Float64}(I,size(Y,2),size(Y,2)),"rrr",k)
    reference=1.0 ./collect(1:length(b)) # external fixed convention; beta[1] is zero
    dot(b,reference)<0 && (b=-b)
    b
end

function pc_worker(planpath,task,out)
    plan=TOML.parsefile(planpath)
    tasks=[(s,r:min(r+plan["shard_size"]-1,s["replicates"])) for s in plan["scenario"] for r in 1:plan["shard_size"]:s["replicates"]]
    s,reps=tasks[task];mkpath(out);path=joinpath(out,"part-$(lpad(task,4,'0')).csv")
    isfile(path) && error("finalized shard exists")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,design,n,ntrain,index,config,C,hfactor,pilot_ratio,adjustment,point,target,a,b,truth,estimate,se,covered,status,reason,index_loss,orientation_flip,missing,extra,effective_n,selected_mean,penalty,seconds")
        p,d=s["p"],s["d"];pairs,labels=pc_targets(s["design"],p)
        for rep in reps
            seed=plan["seed"]+100000*s["id"]+rep
            data=fg_data(seed,s["n"],p,d,s["design"],s["amplitude"])
            perm=randperm(MersenneTwister(seed+39017),p)
            entries=[minmax(invperm(perm)[j],invperm(perm)[k]) for (j,k) in pairs]
            directions=Dict("known"=>data.beta)
            if any(!=("known"),plan["indices"])
                tseed=seed+10000000000
                train=fg_data(tseed,s["ntrain"],p,d,s["design"],s["amplitude"])
                for index in plan["indices"]
                    index=="known" && continue
                    directions[index]=pc_index(train,perm,tseed,s["bank"],index=="estimated_full" ? d : 10)
                end
            end
            for index in plan["indices"]
                b=directions[index];loss=1-dot(b,data.beta)^2;flip=dot(b,data.beta)<0
                missing=count(j->data.beta[j]!=0&&b[j]==0,1:d);extra=count(j->data.beta[j]==0&&b[j]!=0,1:d)
                for cfg in plan["candidate"]
                    haskey(cfg,"design") && cfg["design"]!=s["design"] && continue
                    h=cfg["hfactor"]*s["n"]^(-0.2)
                    adjustment=cfg["adjustment"]
                    started=time_ns()
                    fit=fit_precision(data.X,data.Z*b;points=plan["points"],entries,
                        bandwidth=h,pilot_bandwidth=h*cfg["pilot_ratio"],pilot=:sqrtlasso,
                        refit=true,degree=2,correction=:block,multiplier=cfg["C"],
                        penalty_rule=:normal,penalty_probability=0.1,penalty_tests=1,
                        residual_adjustment=adjustment=="df" ? :df : :none,
                        studentization=adjustment=="hc3" ? :hc3 : :hc0,
                        inverse_bias=adjustment=="inverse" ? :secondorder : :none)
                    seconds=(time_ns()-started)/1e9
                    for (g,t) in enumerate(plan["points"]),(e,(j,k)) in enumerate(pairs)
                        truth=fg_curve(s["design"],p,s["amplitude"],t).O[j,k]
                        value,se=fit.estimate[g,e],fit.se[g,e]
                        reason=replace(fit.failure[g,e],','=>';','\n'=>' ')
                        ok=isfinite(value)&&isfinite(se)&&se>0&&isempty(reason)
                        covered=ok&&abs(value-truth)<=1.959963984540054se
                        diag=filter(q->q.point==t&&q.entry==entries[e],fit.settings.nuisance_diagnostics)
                        neff=isempty(diag) ? NaN : diag[1].effective_n
                        selected=isempty(diag) ? NaN : sum(length,diag[1].supports)/length(diag[1].supports)
                        penalty=isempty(diag) ? NaN : diag[1].penalty
                        println(io,join((s["id"],rep,s["design"],s["n"],s["ntrain"],index,cfg["id"],cfg["C"],cfg["hfactor"],cfg["pilot_ratio"],adjustment,t,labels[e],entries[e][1],entries[e][2],truth,value,se,Int(covered),ok ? "ok" : "failed",reason,loss,Int(flip),missing,extra,neff,selected,penalty,seconds),','))
                    end
                end
            end
            flush(io)
        end
    end
    mv(path*".tmp",path)
end

if abspath(PROGRAM_FILE)==@__FILE__
    pc_worker(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
