"""Replay selected rows only. Primary direction estimates are not recomputed."""
include("multiprobe_common.jl")

function main(planfile,task,out)
    plan=TOML.parsefile(planfile)
    tasks=[(s,first:min(first+plan["shard_size"]-1,s["replicates"]))
           for s in plan["scenario"] for first in 1:plan["shard_size"]:s["replicates"]]
    s,reps=tasks[task]
    configs=haskey(plan,"candidate") ?
        [(id=c["id"],bank=c["bank"],L=c["L"],method=c["method"],k=c["k"]) for c in plan["candidate"]] : mp_configs()
    configs=filter(c->c.method=="screen",configs)
    banks=unique((c.bank,c.L) for c in configs)
    mkpath(out);path=joinpath(out,"part-$(lpad(task,4,'0')).csv")
    isfile(path) && error("finalized support shard exists")
    open(path*".tmp","w") do io
        println(io,"scenario,replicate,id,k,missing,extra,qr_checked,qr_error,direction_leak")
        for rep in reps
            seed=plan["seed"]+100000*s["id"]+rep
            data=mp_data(MersenneTwister(seed),s["n"],s["p"],s["design"];
                signal=get(s,"signal",1.0),beta_kind=get(s,"beta_kind","equal"),even=get(s,"even",false))
            Z=data.Z;Zc=Z.-sum(Z;dims=1)/size(Z,1)
            factor=cholesky(Symmetric(Zc'Zc))
            for (bank,L) in banks
                Y=mp_responses(data.X,mp_bank(s["p"],bank,L,seed+47017))
                mu,scale=CovarianceSelection.SingleIndexGraph._probe_scale(Y)
                W=(Y.-mu)./scale;Wc=W.-sum(W;dims=1)/size(W,1)
                B=factor\(Zc'Wc)
                scores=vec(sum(abs2,B;dims=2));ranking=sortperm(scores;rev=true)
                cs=filter(c->c.bank==bank&&c.L==L,configs)
                ambiguous=any(c->scores[ranking[c.k]]-scores[ranking[c.k+1]]<
                    1e-8*maximum(scores),cs)
                checked=rep<=2 || ambiguous
                err=NaN
                if checked
                    direct=Zc\Wc;err=norm(B-direct)/max(norm(direct),eps(Float64))
                    @assert err<1e-7
                    direct_rank=sortperm(vec(sum(abs2,direct;dims=2));rev=true)
                    if !ambiguous
                        @assert all(sort(ranking[1:c.k])==sort(direct_rank[1:c.k]) for c in cs)
                    end
                    # Direct QR is authoritative in the unlikely near-tie case.
                    ranking=direct_rank
                end
                for c in cs
                    selected=sort(ranking[1:c.k])
                    missing=count(j->data.beta[j]!=0&&!(j in selected),1:20)
                    extra=count(j->data.beta[j]==0&&j in selected,1:20)
                    @assert extra==c.k-4+missing
                    leak=NaN
                    if checked
                        R=Zc[:,selected]\Wc;P=zeros(20,size(W,2));P[selected,:]=R
                        old=CovarianceSelection.SingleIndexGraph._probe_leading(P*P')
                        projected=zeros(20);projected[selected]=old[selected];projected./=norm(projected)
                        leak=norm(old-projected)
                        @assert leak<1e-8
                    end
                    println(io,join((s["id"],rep,c.id,c.k,missing,extra,checked,err,leak),','))
                end
            end
        end
    end
    mv(path*".tmp",path)
end
main(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
