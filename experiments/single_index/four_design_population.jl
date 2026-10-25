using Test
include("four_design_common.jl")
function main(out)
mkpath(out)
xs=collect(range(-7.,7.;length=1201));weights=exp.(-xs.^2/2);weights[[1,end]].*=0.5;weights./=sum(weights)
checks=0
open(joinpath(out,"population.csv"),"w") do io
    println(io,"design,p,amplitude,min_sigma_eigenvalue,max_sigma_eigenvalue,max_condition,first_derivative_fro,scaled_first_derivative_fro,second_stein_fro,mean_gaussian_information,trace_probe_first_signal,oracle_probe_first_signal,stein_identity_error,finite_difference_error")
    for design in ("correlation","opposed_scale","fixed_diagonal","noncommuting"),p in (20,40),a in (design=="correlation" ? (0.15,0.35,0.55) : (0.2,0.4,0.6))
        graphs=fg_graphs(p);bar=zeros(p,p);first=zeros(p,p);stein=zeros(p,p);second=zeros(p,p);info=0.;lo=Inf;hi=0.;cond=0.;fd=0.
        for (t,w) in zip(vcat(xs,[-Inf,Inf]),vcat(weights,[0.,0.]))
            C=fg_curve(design,p,a,t,graphs)
            ev=eigvals(Symmetric(C.S));lo=min(lo,minimum(ev));hi=max(hi,maximum(ev));cond=max(cond,maximum(ev)/minimum(ev))
            bar .+=w.*C.S;first .+=w.*C.D
            if w>0
                stein .+=w*t.*C.S;second .+=w*(t*t-1).*C.S
            end
            info+=w*tr(C.O*C.D*C.O*C.D)/2
        end
        for t in (-1.3,0.,0.8)
            h=1e-5;D=(fg_curve(design,p,a,t+h,graphs).S-fg_curve(design,p,a,t-h,graphs).S)/(2h)
            fd=max(fd,norm(D-fg_curve(design,p,a,t,graphs).D)/max(1.,norm(D)))
        end
        E=eigen(Symmetric(bar));R=E.vectors*Diagonal(1 ./sqrt.(E.values))*E.vectors'
        oracle=-fg_curve(design,p,a,0.,graphs).DO
        err=norm(first-stein)/max(1.,norm(first))
        @test lo>=0.2-1e-10;@test err<1e-7;@test fd<1e-7
        @test norm(first)>1e-3;@test info>0;@test tr(oracle*first)>0
        design=="noncommuting" && @test norm(graphs.A*graphs.B-graphs.B*graphs.A)>0.1
        design in ("correlation","opposed_scale") && @test abs(tr(first))<1e-8
        checks+=6+(design=="noncommuting")+(design in ("correlation","opposed_scale"))
        println(io,join((design,p,a,lo,hi,cond,norm(first),norm(R*first*R),norm(second),info,tr(first),tr(oracle*first),err,fd),','))
    end
end
println("Population diagnostics passed: $checks checks, 24 designs, quadrature on 1201 nodes, analytic endpoint checks.")

end
main(ARGS[1])
