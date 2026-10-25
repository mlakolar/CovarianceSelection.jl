"""
    probe_index_moments(Y, Z; order=1)

Unbiased sample covariance between scalar response probes (columns of `Y`)
and Gaussian Stein features. `Z` must have independent standard Gaussian rows.
Order one returns a d-by-L matrix. Order two returns d-by-d-by-L symmetric
matrices, using the features `Z_j*Z_k - (j==k)`. Constant terms disappear after
centering the responses. This does not estimate the covariate score.
"""
function probe_index_moments(Y::AbstractMatrix,Z::AbstractMatrix;order::Integer=1)
    n,d=size(Z)
    n>=3 && size(Y,1)==n && size(Y,2)>0 || throw(ArgumentError("invalid probe dimensions"))
    all(isfinite,Y) && all(isfinite,Z) || throw(ArgumentError("nonfinite probe data"))
    C=Float64.(Y) .- sum(Y;dims=1)/n
    if order==1
        return Matrix(Z'C/(n-1))
    end
    order==2 || throw(ArgumentError("Stein order must be one or two"))
    pairs=[(j,k) for k in 1:d for j in 1:k]
    F=hcat([Z[:,j].*Z[:,k] for (j,k) in pairs]...)
    V=F'C/(n-1)
    H=zeros(d,d,size(Y,2))
    for (r,(j,k)) in enumerate(pairs)
        H[j,k,:]=V[r,:]
        H[k,j,:]=V[r,:]
    end
    H
end

function _probe_scale(Y)
    n=size(Y,1)
    mu=sum(Y;dims=1)/n
    scale=sqrt.(sum(abs2,Y .- mu;dims=1)/(n-1))
    all(x->isfinite(x)&&x>sqrt(eps(Float64)),scale) ||
        throw(DomainError(scale,"constant or invalid probe"))
    mu,scale
end

function _probe_leading(G;reference=nothing)
    E=eigen(Symmetric((G+G')/2))
    E.values[end]>eps(Float64)*max(1.0,opnorm(G,Inf)) ||
        throw(DomainError(E.values[end],"no positive empirical signal eigenvalue"))
    _orient_index_direction(E.vectors[:,end];reference)
end

"""
    estimate_probe_index(Y, Z; method=:ols, support_size=size(Z,2),
                         split=nothing, standardize=true, reference=nothing)

Estimate a single Gaussian index from fixed scalar response probes without
estimating their conditional means. Methods are `:ols`, `:moment`, `:screen`,
`:split`, `:second`, and `:second_split`. OLS includes an intercept. Screening
ranks rows of the full OLS slope matrix by Euclidean norm, retains a specified
number, then refits OLS before taking a leading singular direction.

Split methods cross-multiply independent-half moments and take the eigenvector
of the largest positive eigenvalue. Order-two aggregation sums products across
probes. By default standardized split methods reserve the first fifth of rows
for independent probe calibration, then split the remaining rows in half.
A supplied split requires `standardize=false` and must partition all rows into
two sets of at least three. Fixed or independently learned scales preserve the
conditional cross-moment identity; calibrating on a moment half would not.
No exact-support or first-coordinate requirement is imposed. Returned directions
have unit norm. These procedures have no automatic inference guarantee.
"""
function estimate_probe_index(Y::AbstractMatrix,Z::AbstractMatrix;
        method::Symbol=:ols,support_size::Integer=size(Z,2),split=nothing,
        standardize::Bool=true,reference=nothing)
    n,d=size(Z)
    size(Y,1)==n && n>d+2 || throw(ArgumentError("need matching rows and n>d+2"))
    all(isfinite,Y) && all(isfinite,Z) || throw(ArgumentError("nonfinite probe data"))
    1<=support_size<=d || throw(ArgumentError("invalid support size"))
    splitmethod=method in (:split,:second_split)
    if splitmethod
        if standardize
            split===nothing || throw(ArgumentError("explicit split requires external scaling"))
            nc=max(3,n÷5)
            calibration=collect(1:nc)
            mid=nc+(n-nc)÷2
            a,b=collect(nc+1:mid),collect(mid+1:n)
        else
            a,b=split===nothing ? (collect(1:n÷2),collect(n÷2+1:n)) : split
            length(a)>=3 && length(b)>=3 && sort(vcat(a,b))==collect(1:n) ||
                throw(ArgumentError("split must partition observations"))
        end
    end
    mu,scale=standardize ? _probe_scale(splitmethod ? Y[calibration,:] : Y) :
        (zeros(1,size(Y,2)),ones(1,size(Y,2)))
    W=(Y .- mu)./scale
    selected=collect(1:d)
    if method in (:ols,:screen)
        Zc=Z .- sum(Z;dims=1)/n
        Wc=W .- sum(W;dims=1)/n
        B=Zc\Wc
        if method==:screen
            selected=sort(sortperm(vec(sum(abs2,B;dims=2));rev=true)[1:support_size])
            R=Zc[:,selected]\Wc
            B=zeros(d,size(W,2))
            B[selected,:]=R
        end
        G=B*B'
    elseif method==:moment
        M=probe_index_moments(W,Z)
        G=M*M'
    elseif method==:split
        M1=probe_index_moments(W[a,:],Z[a,:])
        M2=probe_index_moments(W[b,:],Z[b,:])
        G=(M1*M2'+M2*M1')/2
    elseif method==:second
        H=probe_index_moments(W,Z;order=2)
        G=sum(H[:,:,l]*H[:,:,l] for l in axes(H,3))
    elseif method==:second_split
        H1=probe_index_moments(W[a,:],Z[a,:];order=2)
        H2=probe_index_moments(W[b,:],Z[b,:];order=2)
        G=sum((H1[:,:,l]*H2[:,:,l]+H2[:,:,l]*H1[:,:,l])/2 for l in axes(H1,3))
    else
        throw(ArgumentError("unknown probe index method"))
    end
    direction=_probe_leading(G;reference)
    (;direction,selected,method,support_size,center=vec(mu),scale=vec(scale))
end
