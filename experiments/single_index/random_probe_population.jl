# Average single-probe population signal, without fitting or selecting using truth.
include("four_design_common.jl")
using Test, Statistics
function normal_rule(m)
    e=eigen(SymTridiagonal(zeros(m),sqrt.(collect(1.0:m-1))))
    e.values,vec(e.vectors[1,:].^2)
end
function population(design,p,a,m)
    ts,ws=normal_rule(m);bar=zeros(p,p);D=zeros(p,p);V=zeros(p,p)
    curves=Matrix{Float64}[]
    for (t,w) in zip(ts,ws)
        c=fg_curve(design,p,a,t);push!(curves,c.S)
        bar .+= w.*c.S;D .+=w.*c.D
        V .+=w.*(diag(c.S)*diag(c.S)' .+2 .*c.S.^2)
    end
    V .-=bar.^2
    (;ts,ws,bar,D,V,curves)
end
function random_metrics(P,A,B)
    L=size(A,2);v=zeros(L);mu=zeros(L)
    for (S,w) in zip(P.curves,P.ws)
        SA=S*A;SB=S*B
        aa=vec(sum(A.*SA;dims=1));bb=vec(sum(B.*SB;dims=1));ab=vec(sum(A.*SB;dims=1))
        mu .+=w.*ab;v .+=w.*(aa.*bb.+2 .*ab.^2)
    end
    v .-=mu.^2
    c=vec(sum(A.*(P.D*B);dims=1))
    c.^2,v
end
function summarize(io,design,p,kind,signal2,variance;exact=false)
    r=signal2./variance;information=r./(1 .-r)
    @test minimum(variance)>0;@test minimum(r)>-1e-12;@test maximum(r)<1
    n=length(r);mu=sum(r)/n;mcse=exact ? 0. : std(r)/sqrt(n)
    q=sort(r)
    println(io,join((design,p,kind,n,mu,mcse,q[clamp(ceil(Int,.1n),1,n)],q[clamp(ceil(Int,.5n),1,n)],q[clamp(ceil(Int,.9n),1,n)],sum(information)/n,sum(signal2)/n,count(<(1e-16),signal2)/n,exact),','))
end
function main(out)
    BLAS.set_num_threads(1);mkpath(out);rng=MersenneTwister(65193717);L=20000;checks=0
    open(joinpath(out,"summary.csv"),"w") do io
        println(io,"design,p,probe,count,mean_squared_correlation,mcse,q10,median,q90,mean_directional_information,mean_raw_squared_signal,near_zero_signal_fraction,enumerated")
        for p in (20,40),design in ("opposed_scale","noncommuting","correlation","fixed_diagonal")
            a=design=="correlation" ? .55 : .6
            P=population(design,p,a,121);Q=population(design,p,a,201)
            @test norm(P.D-Q.D)/max(1.,norm(Q.D))<1e-7
            @test norm(P.V-Q.V)/norm(Q.V)<1e-7
            @test abs(sum(P.ws)-1)<1e-12
            @test abs(sum(P.ws.*P.ts.^2)-1)<1e-12
            checks+=4
            js=collect(1:p);idx=[CartesianIndex(j,k) for j in 1:p for k in j+1:p]
            summarize(io,design,p,"coordinate_square",diag(P.D).^2,diag(P.V);exact=true)
            summarize(io,design,p,"coordinate_pair",P.D[idx].^2,P.V[idx];exact=true)
            A=randn(rng,p,L);B=randn(rng,p,L)
            A ./=sqrt.(sum(abs2,A;dims=1));B ./=sqrt.(sum(abs2,B;dims=1))
            signal2,v=random_metrics(P,A,B)
            theory=sum(abs2,P.D)/p^2
            @test abs(sum(signal2)/L-theory)<6std(signal2)/sqrt(L)+1e-12
            checks+=1
            summarize(io,design,p,"two_random_directions",signal2,v)
            signal2,v=random_metrics(P,A,A)
            theory=(tr(P.D)^2+2sum(abs2,P.D))/(p*(p+2))
            @test abs(sum(signal2)/L-theory)<6std(signal2)/sqrt(L)+1e-12
            checks+=1
            summarize(io,design,p,"random_direction_square",signal2,v)
            println("completed ",design," p=",p)
        end
    end
    # Cross-check the general variance formula against Sigma=I.
    A=randn(rng,5,300);B=randn(rng,5,300)
    A ./=sqrt.(sum(abs2,A;dims=1));B ./=sqrt.(sum(abs2,B;dims=1))
    P=(curves=[Matrix{Float64}(I,5,5)],ws=[1.],D=zeros(5,5))
    c,v=random_metrics(P,A,B)
    @test maximum(abs.(v .-(1 .+vec(sum(A.*B;dims=1)).^2)))<1e-12
    checks+=1
    open(joinpath(out,"receipt.txt"),"w") do io
        println(io,"passed $checks checks; 32 summary rows; 20000 random probes per type/configuration; exact enumeration of coordinate probes; GH121 versus GH201 convergence; no index fits or inference")
    end
end
main(ARGS[1])
