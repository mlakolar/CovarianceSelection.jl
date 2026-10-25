# Small paired diagnosis, not a replacement coverage study. Permuting Gaussian
# covariate labels gives a zero true first coefficient with identical X and t.
include("gaussian_fixed_worker.jl")
using Test
include("../../test/index_orientation.jl")
@testset "OLS orientation after support selection" begin
    z=collect(range(-1.,1.;length=120))
    D=(Z=hcat(cos.(z),z,z.^2),X=hcat(ones(120),z))
    @test gfols(D,[2])≈[0.,1.,0.]
    @test_throws DomainError gfols(D,[2];anchor=1)
end

function orientation_check(out)
    mkpath(out)
    open(joinpath(out,"comparison.csv"),"w") do io
        gfrow(io,["n","zero_first","rep","policy","failure_kind","error2","projective_error2","missed","extra","seed"])
        for (id,n) in enumerate((800,3200)),rep in 1:500
            seed=981000000+id*100000+rep
            c=Dict("n"=>n,"p"=>40,"d"=>20,"graph"=>"fixed_chain","noise"=>"gaussian",
                "design"=>"gaussian","index_sparsity"=>4,"index_equal"=>true)
            original=generate_data(seed,c)
            for zero_first in (false,true)
                perm=zero_first ? [5,1,2,3,4,collect(6:20)...] : collect(1:20)
                D=(X=original.X,Z=original.Z[:,perm],score=original.score[:,perm],gamma=original.gamma[perm])
                idx=nothing
                try
                    idx=estimate_index(D.X,D.score;blocks=1,multiplier=1.25,unshrink=true)
                catch err
                    err isa DomainError || rethrow()
                end
                for policy in ("fixed_first","automatic")
                    failure=idx===nothing ? "empty_support" : ""
                    gamma=nothing; selected=idx===nothing ? Int[] : findall(!iszero,idx.direction)
                    if isempty(failure)
                        if policy=="fixed_first" && idx.direction[1]==0
                            failure="missing_anchor"
                        else
                            try
                                gamma=gfols(D,selected;anchor=policy=="fixed_first" ? 1 : nothing)
                            catch err
                                err isa InterruptException && rethrow()
                                failure="ols: "*sprint(showerror,err)
                            end
                        end
                    end
                    truth=findall(!iszero,D.gamma)
                    e=gamma===nothing ? NaN : sum(abs2,gamma-D.gamma)
                    pe=gamma===nothing ? NaN : min(e,sum(abs2,gamma+D.gamma))
                    gfrow(io,[n,Int(zero_first),rep,policy,failure,e,pe,
                        length(setdiff(truth,selected)),length(setdiff(selected,truth)),seed])
                end
            end
        end
    end
end

orientation_check(ARGS[1])
