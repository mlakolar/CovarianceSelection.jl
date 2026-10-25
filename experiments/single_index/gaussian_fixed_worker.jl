# A staged study: index selection is frozen before second-stage simulations.
using CovarianceSelection, LinearAlgebra, Random, TOML
using CovarianceSelection.SingleIndexGraph
const GFSIG = CovarianceSelection.SingleIndexGraph
include("dgp.jl")
BLAS.set_num_threads(1)

function gfrow(io, items)
    println(io,join([x isa AbstractString ? "\""*replace(x,'"'=>"\"\"",'\n'=>' ')*"\"" : string(x) for x in items],','))
end

function gfindex(D,method,multiplier)
    method == "known" && return D.gamma
    method == "oracle_support_ols" && return gfols(D,findall(!iszero,D.gamma);anchor=1)
    blocks = startswith(method,"mom") ? 5 : 1
    idx = estimate_index(D.X,D.score;blocks,multiplier,anchor=1,
                         unshrink=!endswith(method,"soft"))
    endswith(method,"ols") ? gfols(D,findall(!iszero,idx.direction);anchor=1) : idx.direction
end

# Gaussian-Z support refit: the population best linear coefficient for the
# fixed probe response product is parallel to the index. No intercept penalty;
# include an intercept to avoid finite-sample mean contamination.
function gfols(D,support;anchor=nothing,reference=nothing)
    isempty(support) && error("empty index support")
    U = hcat(ones(size(D.Z,1)),D.Z[:,support])
    b = cholesky(Symmetric(U'U)) \ (U'*(D.X[:,end-1].*D.X[:,end]))
    g = zeros(size(D.Z,2))
    g[support] = b[2:end]
    GFSIG._orient_index_direction(g;anchor,reference)
end

function gfshard(planpath,taskid,outdir)
    plan = TOML.parsefile(planpath)
    tasks = [(c,r:min(r+Int(plan["shard_size"])-1,Int(c["replicates"])))
             for c in plan["scenario"] for r in 1:Int(plan["shard_size"]):Int(c["replicates"])]
    c,reps = tasks[taskid]
    mkpath(outdir)
    path = joinpath(outdir,"shard-$(lpad(taskid,5,'0')).csv")
    isfile(path) && error("refusing to overwrite completed shard")
    indexstage = plan["stage"] == "index"
    open(path*".partial","w") do io
        if indexstage
            gfrow(io,["scenario","rep","method","index_multiplier","error2","missed","extra","exact","contains","seconds","failure","seed"])
        else
            gfrow(io,["scenario","rep","index_method","index_multiplier","C","h_factor","pilot_factor","variant","kind","entry_a","entry_b","point","estimate","se","truth","covered","width","failure","index_error2","index_missed","index_extra","index_exact","nuisance_missed","nuisance_extra","nuisance_size","effective_n","lambda","seconds","seed","final_effective_n","final_local_n"])
        end
        for rep in reps
            seed = Int(plan["seed_base"])+Int(c["id"])*100000+rep
            D = generate_data(seed,c)
            for method in c["index_methods"],im in c["index_multipliers"]
                method in ("known","oracle_support_ols") && im != first(c["index_multipliers"]) && continue
                start = time()
                gamma,failure = nothing,""
                try
                    gamma = gfindex(D,method,Float64(im))
                catch err
                    err isa InterruptException && rethrow()
                    failure = sprint(showerror,err)
                end
                secs = time()-start
                tru = Set(findall(!iszero,D.gamma))
                sel = gamma === nothing ? Set{Int}() : Set(findall(!iszero,gamma))
                missed,extra = length(setdiff(tru,sel)),length(setdiff(sel,tru))
                error2 = gamma === nothing ? NaN : sum(abs2,gamma-D.gamma)
                if indexstage
                    gfrow(io,[c["name"],rep,method,im,error2,missed,extra,
                              Int(isempty(failure)&&sel==tru),Int(isempty(failure)&&missed==0),secs,failure,seed])
                    continue
                end
                ti = gamma === nothing ? D.t : D.Z*gamma
                grid = collect(range(Float64(c["grid_min"]),Float64(c["grid_max"]);length=Int(c["grid_n"])))
                entries = [Tuple(Int.(pair)) for pair in get(c,"entries",[[1,1],[1,2],[1,4]])]
                covariances = Dict(point=>D.truth(point) for point in grid)
                precisions = Dict(point=>inv(Symmetric(covariances[point])) for point in grid)
                truth = [precisions[t][a,b] for t in grid,(a,b) in entries]
                for C in c["constants"],hf in c["h_factors"],pf in c["pilot_factors"],variant in c["variants"]
                    h,hb = Float64(hf)*2.5Int(c["n"])^(-0.2),Float64(pf)*2.5Int(c["n"])^(-0.2)
                    fit,fitfailure = nothing,failure
                    union_diagnostics = Dict()
                    start = time()
                    try
                        isempty(failure) || error(failure)
                        oracle = if variant == "oracle_nuisance"
                            (point,J,A)->covariances[point][A,A]\covariances[point][A,J]
                        elseif variant == "union_refit"
                            (point,J,A)->begin
                                initial = GFSIG._regression(D.X,ti,point,hb,J,A,Float64(C),biweight;
                                    refit=true,penalty_rule=:normal,penalty_probability=Float64(plan["penalty_probability"]))
                                union_diagnostics[(point,Tuple(J))] = (point=point,entry=(first(J),last(J)),targets=J,
                                    supports=[A[findall(!iszero,initial.coefficients[:,j])] for j in eachindex(J)],
                                    effective_n=initial.effective_n,penalty=initial.penalty)
                                GFSIG._residual_df(D.X,ti,point,hb,J,A,initial.coefficients,biweight).coefficients
                            end
                        else
                            nothing
                        end
                        fit = fit_precision(D.X,ti;points=grid,entries,bandwidth=h,pilot_bandwidth=hb,
                            pilot=:sqrtlasso,refit=true,multiplier=Float64(C),penalty_rule=:normal,
                            penalty_probability=Float64(plan["penalty_probability"]),
                            residual_adjustment=variant=="corrected" ? :df : :none,
                            studentization=variant=="corrected" ? :hc3 : :hc0,
                            inverse_bias=variant=="corrected" ? :secondorder : :none,
                            oracle_coefficients=oracle)
                    catch err
                        err isa InterruptException && rethrow()
                        fitfailure = sprint(showerror,err)
                    end
                    seconds = time()-start
                    ci = fit === nothing ? nothing : confint(fit)
                    for (e,(a,b)) in enumerate(entries),(g,point) in enumerate(grid)
                        fail = fit === nothing ? fitfailure : fit.failure[g,e]
                        valid = isempty(fail)
                        diag = fit === nothing ? [] : filter(q->q.point==point && q.entry==(a,b),fit.settings.nuisance_diagnostics)
                        if variant=="union_refit" && fit!==nothing
                            Jkey = a==b ? (a,) : (a,b)
                            diag = haskey(union_diagnostics,(point,Jkey)) ? [union_diagnostics[(point,Jkey)]] : []
                        end
                        nm,ne,ns,en,lambda = NaN,NaN,NaN,NaN,NaN
                        if !isempty(diag)
                            q = only(diag)
                            J = q.targets; A = setdiff(1:Int(c["p"]),J)
                            B = covariances[point][A,A]\covariances[point][A,J]
                            # Report union support: the corrected variant refits
                            # this union; the plain variant refits columns separately.
                            tr = Set(A[findall(vec(any(abs.(B).>1e-10;dims=2)))])
                            ss = Set(vcat(q.supports...))
                            nm,ne,ns = length(setdiff(tr,ss)),length(setdiff(ss,tr)),length(ss)
                            en,lambda = q.effective_n,q.penalty
                        end
                        cover = valid && ci.lower[g,e]<=truth[g,e]<=ci.upper[g,e]
                        gfrow(io,[c["name"],rep,method,im,C,hf,pf,variant,"pointwise",a,b,point,
                            valid ? fit.estimate[g,e] : NaN,valid ? fit.se[g,e] : NaN,truth[g,e],Int(cover),
                            valid ? ci.upper[g,e]-ci.lower[g,e] : NaN,fail,error2,missed,extra,Int(isempty(failure)&&sel==tru),
                            nm,ne,ns,en,lambda,seconds,seed,
                            fit===nothing ? NaN : fit.diagnostics[g].equivalent_n,
                            fit===nothing ? 0 : fit.diagnostics[g].local_n])
                    end
                    if Int(c["draws"])>=100
                        for joint in (false,true)
                            band,fail = nothing,fitfailure
                            try
                                fit===nothing && error(fitfailure)
                                band=confband(fit;draws=Int(c["draws"]),joint,rng=MersenneTwister(seed+800000000))
                            catch err
                                err isa InterruptException && rethrow()
                                fail=sprint(showerror,err)
                            end
                            for (e,(a,b)) in enumerate(entries)
                                covered=band!==nothing && all(band.lower[:,e].<=truth[:,e].<=band.upper[:,e])
                                width=band===nothing ? NaN : sum(band.upper[:,e]-band.lower[:,e])/length(grid)
                                gfrow(io,[c["name"],rep,method,im,C,hf,pf,variant,joint ? "family_component" : "curve",a,b,NaN,
                                    NaN,NaN,NaN,Int(covered),width,fail,error2,missed,extra,Int(isempty(failure)&&sel==tru),
                                    NaN,NaN,NaN,NaN,NaN,seconds,seed,NaN,NaN])
                            end
                            if joint
                                covered=band!==nothing && all(band.lower.<=truth.<=band.upper)
                                width=band===nothing ? NaN : sum(band.upper-band.lower)/length(truth)
                                gfrow(io,[c["name"],rep,method,im,C,hf,pf,variant,"family",0,0,NaN,NaN,NaN,NaN,
                                    Int(covered),width,fail,error2,missed,extra,Int(isempty(failure)&&sel==tru),
                                    NaN,NaN,NaN,NaN,NaN,seconds,seed,NaN,NaN])
                            end
                        end
                    end
                end
            end
            flush(io)
            println("stage=",plan["stage"]," scenario=",c["name"]," rep=",rep)
            flush(stdout)
        end
    end
    mv(path*".partial",path)
end

if abspath(PROGRAM_FILE)==@__FILE__
    gfshard(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
