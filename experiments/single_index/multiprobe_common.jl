using LinearAlgebra, Random, TOML
import CovarianceSelection
using CovarianceSelection.SingleIndexGraph

const MP_METHODS=[("ols",20),("moment",20),("screen",6),("screen",10),
                  ("split",20),("second",20),("second_split",20)]

function mp_beta(kind)
    b=zeros(20)
    b[2:5]=kind=="mixed" ? [0.7,-0.5,0.4,-0.3] : ones(4)
    b/norm(b)
end

mp_link(t,even)=even ? tanh(t)^2 : tanh(t)+0.5tanh(t)^2

function mp_data(rng,n,p,design;signal=1.0,beta_kind="equal",even=false)
    beta=mp_beta(beta_kind)
    Z=randn(rng,n,20)
    t=Z*beta
    g=signal.*mp_link.(t,even)
    E=randn(rng,n,p)
    if design=="chain"
        V=[sqrt(2/(p+1))*sin(j*k*pi/(p+1)) for j in 1:p,k in 1:p]
        lam=[2cos(k*pi/(p+1)) for k in 1:p]
        D=1.5 .+0.15g
        edge=0.30 .+0.08g
        X=(E./sqrt.(D .- edge*lam'))*V'
    elseif design=="pairs"
        rho=0.35 .+0.15g
        X=copy(E)
        for j in 1:2:p
            X[:,j+1]=rho.*E[:,j].+sqrt.(1 .-rho.^2).*E[:,j+1]
        end
    else
        error("unknown design")
    end
    (;X,Z,beta)
end

function mp_covariance(t,p,design;signal=1.0,even=false)
    g=signal*mp_link(t,even)
    if design=="chain"
        return inv(SymTridiagonal(fill(1.5+0.15g,p),fill(-(0.30+0.08g),p-1)))
    end
    S=Matrix{Float64}(I,p,p)
    for j in 1:2:p
        S[j,j+1]=S[j+1,j]=0.35+0.15g
    end
    S
end

function mp_bank(p,kind,L,seed)
    rng=MersenneTwister(seed)
    pairs=[(j,k) for k in 1:p for j in 1:k-1]
    shuffle!(rng,pairs)
    q=Vector{Vector{Float64}}(); r=Vector{Vector{Float64}}()
    unit(j)=(v=zeros(p);v[j]=1;v)
    if kind in ("diagonal","mixed","complete")
        for j in 1:p
            push!(q,unit(j));push!(r,unit(j))
        end
    end
    if kind in ("products","mixed","complete","single")
        count=kind=="complete" ? length(pairs) : kind=="single" ? 1 : kind=="mixed" ? L÷2 : L
        for (j,k) in pairs[1:count]
            push!(q,unit(j));push!(r,unit(k))
        end
    end
    if kind in ("projections","mixed")
        count=kind=="mixed" ? L÷2 : L
        for l in 1:count
            v=zeros(p)
            v[randperm(rng,p)[1:4]]=ifelse.(rand(rng,Bool,4),0.5,-0.5)
            push!(q,v);push!(r,v)
        end
    end
    (;Q=hcat(q...),R=hcat(r...),kind,L)
end

function mp_responses(X,bank)
    # Unit-coordinate probes need elementwise products, not dense products
    # against thousands of almost entirely zero projection columns.
    L=size(bank.Q,2)
    simple=bank.kind=="projections" ? 0 : bank.kind=="mixed" ? size(X,2)+bank.L÷2 : L
    Y=Matrix{Float64}(undef,size(X,1),L)
    for l in 1:simple
        j=argmax(view(bank.Q,:,l));k=argmax(view(bank.R,:,l))
        Y[:,l]=X[:,j].*X[:,k]
    end
    if simple<L
        Y[:,simple+1:end]=(X*bank.Q[:,simple+1:end]).*(X*bank.R[:,simple+1:end])
    end
    Y
end

function mp_configs()
    out=NamedTuple[]
    for (bank,L) in vcat([("single",1),("diagonal",0),("complete",0)],
                         [(b,L) for b in ("products","projections","mixed") for L in (10,30,60)])
        methods=bank=="single" ? [("ols",20)] : bank=="complete" ? MP_METHODS[1:4] : MP_METHODS
        for (method,k) in methods
            push!(out,(;id="$(bank)-$(L)-$(method)-$(k)",bank,L,method,k))
        end
    end
    out
end

mp_fit(Y,Z,c)=estimate_probe_index(Y,Z;method=Symbol(c.method),support_size=c.k)

function mp_metrics(u,beta;selected=eachindex(beta))
    cosine=clamp(abs(dot(u,beta)),0.0,1.0)
    (;loss=1-cosine^2,prediction=2-2cosine,angle=acos(cosine)*180/pi,
      missing=count(j->beta[j]!=0&&!(j in selected),eachindex(beta)),
      extra=count(j->beta[j]==0&&j in selected,eachindex(beta)))
end

function mp_validate(X,Z,configs,seed;folds=3)
    n=size(Z,1); rng=MersenneTwister(seed)
    perm=randperm(rng,n); labels=zeros(Int,n); labels[perm]=mod.(0:n-1,folds).+1
    score=zeros(length(configs)); bad=falses(length(configs))
    # The SAME rich reference dictionary is used for every candidate. Candidate
    # dictionary size cannot inflate the criterion by increasing its dimension.
    reference=mp_bank(size(X,2),"complete",0,seed+917)
    Yref=mp_responses(X,reference)
    for fold in 1:folds
        tr=findall(!=(fold),labels); val=findall(==(fold),labels)
        half=length(val)÷2; a=val[1:half]; b=val[half+1:end]
        mu,scale=CovarianceSelection.SingleIndexGraph._probe_scale(Yref[tr,:])
        V=(Yref .-mu)./scale
        M1=probe_index_moments(V[a,:],Z[a,:]); M2=probe_index_moments(V[b,:],Z[b,:])
        cache=Dict{Tuple{String,Int},Matrix{Float64}}()
        for (i,c) in enumerate(configs)
            try
                Y=get!(cache,(c.bank,c.L)) do
                    mp_responses(X,mp_bank(size(X,2),c.bank,c.L,seed-11))
                end
                u=mp_fit(Y[tr,:],Z[tr,:],c).direction
                score[i]+=dot(M1'u,M2'u)/folds
            catch err
                err isa InterruptException && rethrow()
                bad[i]=true
            end
        end
    end
    score[bad].=-Inf
    maximum(score)>-Inf || error("all validation candidates failed")
    argmax(score),score
end
