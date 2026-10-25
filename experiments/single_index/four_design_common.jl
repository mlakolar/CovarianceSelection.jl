# Four approved Gaussian covariance families. Estimators reuse frozen aggregation code.
include("aggregation_common.jl")

function fg_graphs(p;seed=403519)
    A=Matrix(SymTridiagonal(zeros(p),fill(0.5,p-1)))
    rng=MersenneTwister(seed+p)
    B=zeros(p,p)
    for attempt in 1:10000
        fill!(B,0.);v=randperm(rng,p)
        all(abs(v[j]-v[j+1])>1 for j in 1:2:p) || continue
        for j in 1:2:p;B[v[j],v[j+1]]=B[v[j+1],v[j]]=1.;end
        norm(A*B-B*A)>0.1 && return (;A,B)
    end
    error("unable to construct disjoint noncommuting matching")
end

function fg_curve(design,p,a,t,graphs=fg_graphs(p))
    u=tanh(t);du=1-u*u
    S=zeros(p,p);D=zeros(p,p)
    if design=="correlation"
        r=0.15+a*u;dr=a*du
        for j in 1:2:p
            S[j,j]=S[j+1,j+1]=1.;S[j,j+1]=S[j+1,j]=r
            D[j,j+1]=D[j+1,j]=dr
        end
        O=inv(Symmetric(S));DO=-O*D*O
    elseif design=="opposed_scale"
        for (b,j) in enumerate(1:2:p)
            sign=isodd(b) ? 1. : -1.;v=exp(sign*a*u);dv=sign*a*du*v
            S[j:j+1,j:j+1]=v.*[1. 0.25;0.25 1.]
            D[j:j+1,j:j+1]=dv.*[1. 0.25;0.25 1.]
        end
        O=inv(Symmetric(S));DO=-O*D*O
    else
        O=Matrix{Float64}(I,p,p)
        if design=="fixed_diagonal"
            O .+=(0.1+a*u).*graphs.A;DO=a*du.*graphs.A
        elseif design=="noncommuting"
            O .+=0.2 .*graphs.A+a*u.*graphs.B;DO=a*du.*graphs.B
        else
            error("unknown design")
        end
        S=inv(Symmetric(O));D=-S*DO*S
    end
    (;S=Matrix(S),D=Matrix(D),O=Matrix(O),DO=Matrix(DO))
end

function fg_beta(d)
    b=zeros(d)
    b[d==4 ? (1:4) : (2:5)]=[1.,1.,-1.,-1.]/2
    b
end

function fg_data(seed,n,p,d,design,a)
    rng=MersenneTwister(seed);beta=fg_beta(d);Z=randn(rng,n,d);t=Z*beta
    E=randn(rng,n,p);X=similar(E);graphs=fg_graphs(p)
    if design in ("correlation","opposed_scale")
        for (b,j) in enumerate(1:2:p)
            rho=design=="correlation" ? 0.15 .+a.*tanh.(t) : fill(0.25,n)
            scale=design=="opposed_scale" ? exp.((isodd(b) ? 1. : -1.)*a.*tanh.(t)./2) : ones(n)
            X[:,j]=scale.*E[:,j]
            X[:,j+1]=scale.*(rho.*E[:,j]+sqrt.(1 .-rho.^2).*E[:,j+1])
        end
    elseif design=="fixed_diagonal"
        eig=eigen(Symmetric(graphs.A));b=0.1 .+a.*tanh.(t)
        X=(E./sqrt.(1 .+b*eig.values'))*eig.vectors'
    else
        for i in 1:n
            O=Matrix{Float64}(I,p,p)+0.2graphs.A+a*tanh(t[i])*graphs.B
            L=cholesky(Symmetric(O)).L
            X[i,:]=L'\view(E,i,:)
        end
    end
    # A fresh independent permutation hides the response graph from the estimators.
    perm=randperm(MersenneTwister(seed+39017),p)
    X=X[:,perm]
    score=-fg_curve(design,p,a,0.,graphs).DO[perm,perm]
    (;X,Z,beta,score)
end

function fg_metrics(u,beta,Zeval,n)
    all(isfinite,u) && abs(norm(u)-1)<1e-7 || error("invalid direction")
    dot(u,beta)<0 && (u=-u)
    c=clamp(dot(u,beta),0.,1.);loss=1-c*c;prediction=2-2c
    delta=norm(u-beta);t=Zeval*beta;that=Zeval*u
    distortions=Float64[]
    for mult in (0.75,1.,1.5)
        h=mult*n^(-0.2)
        w=exp.(-0.5 .*(t./h).^2);w./=sum(w)
        v=exp.(-0.5 .*(that./h).^2);v./=sum(v)
        push!(distortions,sum(abs,w-v))
    end
    selected=findall(!iszero,u)
    (;loss,prediction,delta,ratio=delta/n^(-0.2),small=distortions[1],l1=distortions[2],large=distortions[3],
      missing=count(j->beta[j]!=0&&!(j in selected),eachindex(beta)),extra=count(j->beta[j]==0&&j in selected,eachindex(beta)))
end
