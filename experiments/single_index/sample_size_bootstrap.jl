# Pairs bootstrap repeats the complete feasible procedure. Failed draws are
# retained. A bootstrap interval is unavailable unless every requested draw
# produces a finite pivot for that entry and location.
include("gaussian_fixed_worker.jl")

function ssfit(D,ti,grid,entries,h,hb,C)
    coefficients = (point,J,A)->begin
        initial=GFSIG._regression(D.X,ti,point,hb,J,A,C,biweight;
            refit=true,penalty_rule=:normal,penalty_probability=.1)
        GFSIG._residual_df(D.X,ti,point,hb,J,A,initial.coefficients,biweight).coefficients
    end
    fit_precision(D.X,ti;points=grid,entries,bandwidth=h,pilot_bandwidth=hb,
        pilot=:sqrtlasso,refit=true,multiplier=C,penalty_rule=:normal,
        penalty_probability=.1,oracle_coefficients=coefficients)
end

function ssquantile(x,p)
    y=sort(x); z=1+(length(y)-1)*p; lo=floor(Int,z); hi=ceil(Int,z)
    (hi==lo) ? y[lo] : y[lo]+(z-lo)*(y[hi]-y[lo])
end

function ssbootstrap(planpath,taskid,outdir)
    plan=TOML.parsefile(planpath)
    tasks=[(c,r:min(r+Int(plan["shard_size"])-1,Int(c["replicates"])))
        for c in plan["scenario"] for r in 1:Int(plan["shard_size"]):Int(c["replicates"])]
    c,reps=tasks[taskid]; mkpath(outdir)
    path=joinpath(outdir,"shard-$(lpad(taskid,5,'0')).csv")
    isfile(path) && error("refusing to overwrite completed bootstrap shard")
    grid=collect(range(Float64(c["grid_min"]),Float64(c["grid_max"]);length=Int(c["grid_n"])))
    entries=[Tuple(Int.(pair)) for pair in c["entries"]]
    n=Int(c["n"]); B=Int(plan["bootstrap_draws"])
    h=2.5n^(-.2); hb=2h; C=1.1
    open(path*".partial","w") do io
        gfrow(io,["scenario","rep","index_method","interval_method","entry_a","entry_b","point",
            "estimate","se","truth","covered","lower","upper","width","failure",
            "bootstrap_failures","bootstrap_draws","index_error2","seed"])
        for rep in reps
            seed=Int(plan["seed_base"])+Int(c["id"])*100000+rep
            D=generate_data(seed,c)
            truth=[inv(Symmetric(D.truth(t)))[a,b] for t in grid,(a,b) in entries]
            for method in c["index_methods"]
                gamma,original,basefail=nothing,nothing,""
                try
                    gamma=gfindex(D,method,1.25)
                    original=ssfit(D,D.Z*gamma,grid,entries,h,hb,C)
                catch err
                    err isa InterruptException && rethrow()
                    basefail=sprint(showerror,err)
                end
                pivots=Dict(mode=>fill(NaN,B,length(grid),length(entries)) for mode in ("fixed_index","full_pipeline"))
                if original!==nothing
                    rng=MersenneTwister(seed+700000000)
                    for b in 1:B
                        ids=rand(rng,1:n,n)
                        Db=(X=D.X[ids,:],Z=D.Z[ids,:],score=D.score[ids,:],gamma=D.gamma)
                        for mode in ("fixed_index","full_pipeline")
                            try
                                gstar=mode=="fixed_index" ? gamma : gfindex(Db,method,1.25)
                                fit=ssfit(Db,Db.Z*gstar,grid,entries,h,hb,C)
                                for e in eachindex(entries),g in eachindex(grid)
                                    if isempty(original.failure[g,e]) && isempty(fit.failure[g,e])
                                        pivots[mode][b,g,e]=(fit.estimate[g,e]-original.estimate[g,e])/fit.se[g,e]
                                    end
                                end
                            catch err
                                err isa InterruptException && rethrow()
                            end
                        end
                    end
                end
                for mode in ("wald","fixed_index","full_pipeline"),e in eachindex(entries),g in eachindex(grid)
                    fail=original===nothing ? basefail : original.failure[g,e]
                    theta=original===nothing ? NaN : original.estimate[g,e]
                    se=original===nothing ? NaN : original.se[g,e]
                    missing=mode=="wald" ? 0 : count(!isfinite,view(pivots[mode],:,g,e))
                    lower,upper=NaN,NaN
                    if isempty(fail)
                        if mode=="wald"
                            lower,upper=theta-1.959963984540054se,theta+1.959963984540054se
                        elseif missing==0
                            ts=collect(view(pivots[mode],:,g,e))
                            lower,upper=theta-ssquantile(ts,.975)*se,theta-ssquantile(ts,.025)*se
                        else
                            fail="bootstrap: $missing/$B invalid draws"
                        end
                    end
                    a,b=entries[e]
                    gfrow(io,[c["name"],rep,method,mode,a,b,grid[g],theta,se,truth[g,e],
                        Int(isempty(fail)&&lower<=truth[g,e]<=upper),lower,upper,upper-lower,fail,
                        missing,B,gamma===nothing ? NaN : sum(abs2,gamma-D.gamma),seed])
                end
            end
            flush(io);println("bootstrap scenario=",c["name"]," rep=",rep);flush(stdout)
        end
    end
    mv(path*".partial",path)
end

if abspath(PROGRAM_FILE)==@__FILE__
    ssbootstrap(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
