# Fixed-rule aggregation experiments. No link fitting and no truth-based tuning.
include("multiprobe_common.jl")

function ag_direction(G)
    u=eigen(Symmetric((G+G')/2)).vectors[:,end]
    u/norm(u)
end

function ag_prepare(Y,Z)
    n=size(Z,1)
    Zc=Z .-sum(Z;dims=1)/n
    mu,scale=CovarianceSelection.SingleIndexGraph._probe_scale(Y)
    V=(Y .-mu)./scale
    G=Symmetric(Zc'Zc/(n-1))
    root=cholesky(G).L
    M=Zc'V/(n-1)
    # Z covariance = root*root'. Transform to identity, then back-transform.
    N=root\M
    B=root'\N
    (;Zc,V,M,N,B,root)
end

function ag_fit_matrix(P,W,method,k)
    d=size(P.Zc,2)
    if k==d
        A=method=="ols" ? P.B : method=="moment" ? P.M : P.N
        u=ag_direction(A*W*A')
        return method=="rrr" ? normalize(P.root'\u) : u
    end
    selected=sort(sortperm(vec(sum(abs2,P.B;dims=2));rev=true)[1:k])
    # Refit on the same selected rows; preserve probe scaling and weight matrix.
    Z=P.Zc[:,selected]
    G=cholesky(Symmetric(Z'Z/(size(Z,1)-1))).L
    M=Z'P.V/(size(Z,1)-1)
    A=method=="ols" ? G'\(G\M) : method=="moment" ? M : G\M
    v=ag_direction(A*W*A')
    method=="rrr" && (v=normalize(G'\v))
    u=zeros(d);u[selected]=v
    u
end

function ag_weights(P)
    L=size(P.V,2)
    C=Symmetric(P.V'P.V/(size(P.V,1)-1))
    E=eigen(C)
    values=max.(E.values,0.0)
    [(rho,E.vectors*Diagonal(1 ./((1-rho).*values .+rho))*E.vectors')
       for rho in (0.05,0.2,0.5)]
end

function ag_learned(Y,Z;rho=0.2)
    n,d=size(Z); mid=n÷2
    halves=(collect(1:mid),collect(mid+1:n))
    moments=Vector{Vector{Float64}}();teachers=Vector{Vector{Float64}}()
    for f in 1:2
        tr,te=halves[f],halves[3-f]
        P=ag_prepare(Y[tr,:],Z[tr,:])
        L=size(Y,2)
        C=Symmetric(P.V'P.V/(length(tr)-1))
        teacher=ag_fit_matrix(P,Matrix{Float64}(I,L,L),"rrr",d)
        # Estimate the shared loading and stabilize redundant probe combinations.
        loading=P.M'teacher
        w=cholesky(Symmetric((1-rho).*C+rho*I))\loading
        w ./=sqrt(max(dot(w,C*w),eps(Float64)))
        mu,scale=CovarianceSelection.SingleIndexGraph._probe_scale(Y[tr,:])
        response=((Y[te,:] .-mu)./scale)*w
        Zc=Z[te,:] .-sum(Z[te,:];dims=1)/length(te)
        response .-=sum(response)/length(response)
        m=Zc'response/(length(te)-1)
        push!(moments,m);push!(teachers,teacher)
    end
    # Fix relative signs using training directions, never the true direction.
    dot(teachers[1],teachers[2])<0 && (moments[2].*=-1)
    normalize(moments[1]+moments[2])
end

function ag_oracle_probe(X,design,signal)
    p=size(X,2)
    if design=="chain"
        return signal.*(-0.15.*vec(sum(abs2,X;dims=2)) .+
                         0.16.*vec(sum(X[:,1:end-1].*X[:,2:end];dims=2)))
    end
    rho=0.35;dr=0.15signal
    diagonal=-2rho*dr/(1-rho^2)^2
    off=dr*(1+rho^2)/(1-rho^2)^2
    diagonal.*vec(sum(abs2,X;dims=2)) .+
       2off.*vec(sum(X[:,1:2:end].*X[:,2:2:end];dims=2))
end

function ag_ols_scalar(y,Z;support=collect(1:size(Z,2)))
    A=Z[:,support] .-sum(Z[:,support];dims=1)/size(Z,1)
    b=A\(y .-sum(y)/length(y));u=zeros(size(Z,2));u[support]=b
    normalize(u)
end
