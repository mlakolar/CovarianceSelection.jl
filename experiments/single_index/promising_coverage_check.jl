include("promising_coverage_worker.jl")
using Test
function checks()
    checks=0
    for family in ("opposed_scale","noncommuting","correlation")
        a=family=="correlation" ? .55 : .6;p=40;seed=238151
        data=fg_data(seed,200,p,20,family,a)
        perm=randperm(MersenneTwister(seed+39017),p)
        pairs,labels=pc_targets(family,p)
        @test length(unique(pairs))==4
        @test dot(fg_beta(20),1.0 ./collect(1:20))>0
        @test data.beta[1]==0
        # Permutation map preserves canonical precision target entries.
        O=fg_curve(family,p,a,0.).O;P=O[perm,perm]
        for (j,k) in pairs
            @test P[invperm(perm)[j],invperm(perm)[k]]==O[j,k];checks+=1
        end
        trainseed=seed+10000000000;train=fg_data(trainseed,200,p,20,family,a)
        b=pc_index(train,perm,trainseed,family=="opposed_scale" ? "mixed" : "complete",10)
        @test norm(b)≈1
        @test count(!iszero,b)==10
        @test dot(b,1.0 ./collect(1:20))>=0
        checks+=6
    end
    println("passed ",checks," design/index checks")
end
checks()
