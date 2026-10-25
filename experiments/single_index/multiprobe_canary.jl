include("multiprobe_common.jl")
using Test
include("../../test/probe_index.jl")

function mp_population(out)
    # Gaussian quadrature integrates expectations under a standard normal index.
    E=eigen(SymTridiagonal(zeros(100),sqrt.(Float64.(1:99))))
    nodes=E.values; weights=E.vectors[1,:].^2
    open(joinpath(out,"population.csv"),"w") do io
        println(io,"design,even,bank,probes,first_norm,second_norm,min_eigenvalue")
        for design in ("chain","pairs"),even in (false,true),kind in ("diagonal","mixed","complete")
            bank=mp_bank(40,kind,60,94811); L=size(bank.Q,2)
            c=zeros(L);b=zeros(L);mu=zeros(L); ey2=zeros(L);mineig=Inf
            for (t,w) in zip(nodes,weights)
                S=mp_covariance(t,40,design;even)
                mineig=min(mineig,eigmin(Symmetric(S)))
                SQ=S*bank.Q; SR=S*bank.R
                m=vec(sum(bank.Q.*SR;dims=1))
                v=vec(sum(bank.Q.*SQ;dims=1)).*vec(sum(bank.R.*SR;dims=1)).+m.^2
                c.+=w*t.*m; b.+=w*(t^2-1).*m; mu.+=w.*m; ey2.+=w.*(v.+m.^2)
            end
            scale=sqrt.(ey2.-mu.^2)
            cn,bn=norm(c./scale),norm(b./scale)
            @test mineig>0.1
            if design=="pairs" && kind=="diagonal"
                @test cn<1e-10 && bn<1e-10
            elseif even
                @test cn<1e-10 && bn>1e-5
            else
                @test cn>1e-5 && bn>1e-5
            end
            println(io,join((design,even,kind,L,cn,bn,mineig),','))
        end
    end
end

function main(out)
    mkpath(out)
    mp_population(out)
    open(joinpath(out,"timings.csv"),"w") do io
        println(io,"design,n,p,id,seconds,loss")
        for (di,(design,n,p)) in enumerate([("chain",200,40),("pairs",400,40),
                                           ("chain",3200,40),("pairs",1600,160)])
            data=mp_data(MersenneTwister(110+di),n,p,design)
            for c in mp_configs()
                start=time_ns()
                Y=mp_responses(data.X,mp_bank(p,c.bank,c.L,170+di))
                u=mp_fit(Y,data.Z,c).direction
                loss=mp_metrics(u,data.beta).loss
                @test isfinite(loss) && 0<=loss<=1
                println(io,join((design,n,p,c.id,(time_ns()-start)/1e9,loss),','))
            end
            c=[(id="mixed-30-ols-20",bank="mixed",L=30,method="ols",k=20),
               (id="mixed-30-screen-6",bank="mixed",L=30,method="screen",k=6)]
            start=time_ns()
            chosen,scores=mp_validate(data.X,data.Z,c,181+di)
            @test 1<=chosen<=2 && all(isfinite,scores)
            println(io,join((design,n,p,"validation2",(time_ns()-start)/1e9,NaN),','))
        end
    end
end
main(ARGS[1])
