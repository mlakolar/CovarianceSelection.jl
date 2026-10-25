# Rank-one projected multivariate MAVE. Probe scaling is fixed across iterations.

function _mave_smoothers(t, anchors, h)
    delta=permutedims(t) .- t[anchors]
    W=exp.(-0.5 .* (delta./h).^2)
    W ./= sum(W;dims=2)
    first=sum(W.*delta;dims=2)
    second=sum(W.*delta.^2;dims=2)
    determinant=second.-first.^2
    minimum(determinant)>eps(Float64)*max(h*h,1.0) ||
        throw(DomainError(minimum(determinant),"unidentified local slope"))
    A=W.*(second.-first.*delta)./determinant
    B=W.*(delta.-first)./determinant
    (;W,A,B)
end

function _mave_system(V,Z,b;bandwidth_factor=1.0,trim=0.05,smoother=nothing)
    n,d=size(Z);L=size(V,2);t=Z*b
    h=bandwidth_factor*n^(-0.2)*sqrt(sum(abs2,t .-sum(t)/n)/(n-1))
    h>0 || throw(DomainError(h,"zero projected bandwidth"))
    order=sortperm(t);cut=floor(Int,trim*n)
    anchors=sort(order[cut+1:n-cut]);m=length(anchors)
    localfit=_mave_smoothers(t,anchors,h)
    A,B=smoother===nothing ? (localfit.A,localfit.B) : smoother(t,anchors,h)
    size(A)==size(B)==(m,n) || throw(DimensionMismatch("invalid smoother matrices"))
    all(isfinite,A) && all(isfinite,B) || error("nonfinite smoother")
    intercept=A*V;slope=B*V;W=localfit.W
    c=vec(sum(abs2,slope;dims=2));Za=Z[anchors,:]
    R=slope*V' .-sum(slope.*intercept;dims=2)
    WR=W.*R;WZ=W*Z
    q=W'c
    cross=Za'*(c.*WZ)
    Q=(Z'*(q.*Z)-cross-cross'+Za'*(c.*Za))/(m*L)
    Q=Matrix(Symmetric((Q+Q')/2))
    r=vec(sum(WR*Z .-Za.*sum(WR;dims=2);dims=1))/(m*L)
    (;Q,r,h,anchors,A,B,W,intercept,slope)
end

function _mave_quadratic(Q,r;support=collect(eachindex(r)),penalty=nothing,
                         tolerance=1e-9,maxiter=10000)
    d=length(r);out=zeros(d);S=collect(support)
    isempty(S) && throw(ArgumentError("empty MAVE support"))
    ridge=1e-10*max(sum(diag(Q))/d,eps(Float64))
    H=Matrix(Q[S,S])+ridge*I;v=r[S]
    if penalty===nothing
        out[S]=cholesky(Symmetric(H))\v
        return (;direction=out,ridge,kkt=0.0,iterations=1)
    end
    length(penalty)==d && all(x->isfinite(x)&&x>=0,penalty) ||
        throw(ArgumentError("invalid coordinate penalties"))
    x=zeros(length(S));lambda=penalty[S];residual=copy(v)
    violation=Inf;iteration=0
    for it in 1:maxiter
        iteration=it
        for j in eachindex(x)
            old=x[j];z=residual[j]+H[j,j]*old
            x[j]=copysign(max(abs(z)-lambda[j],0.0),z)/H[j,j]
            residual .-= H[:,j].*(x[j]-old)
        end
        violation=maximum(iszero(x[j]) ? max(abs(residual[j])-lambda[j],0.0) :
                          abs(residual[j]-lambda[j]*sign(x[j])) for j in eachindex(x))
        violation<=tolerance*max(maximum(abs,v),eps(Float64)) && break
    end
    violation<=tolerance*max(maximum(abs,v),eps(Float64)) || error("MAVE Lasso KKT failure")
    out[S]=x
    (;direction=out,ridge,kkt=violation,iterations=iteration)
end

"""
    refine_probe_index(Y, Z, initial; updates=50, bandwidth_factor=1.0,
                       trim=0.05, support=1:size(Z,2), tolerance=1e-4,
                       standardize=true, smoother=nothing)

Refine a unit index using projected, rank-one multivariate MAVE. Each column
of `Y` is a fixed response probe. Fit local linear probe means and derivatives
at the current index, freeze their Gaussian neighborhoods, and solve the joint
quadratic index update. Normalize the unconstrained update; rescaling its local
slopes gives the same fitted products. The next iteration refits all curves.
The changing-neighborhood objective need not decrease between iterations.

`penalty_method=:lasso` or `:adaptive_lasso` penalizes every quadratic update
using `penalty_level * lambda_max` at that update. Adaptive weights are fixed
from the initializer, with a coefficient floor of 0.05. All allowed coordinates
remain eligible at each update. These options return a penalized direction;
use `select_mave_variables(...; iterative=true)` to refit its final support.

`updates=1` performs exactly one update. All covariates can reenter unless
`support` is explicitly restricted. Trim the bottom and top `trim` fractions
of projected anchor points; all observations remain available as neighbors.
The Gaussian bandwidth is `bandwidth_factor * n^(-1/5) * sd(Z*direction)`.
A supplied `smoother(t, anchors, h)` returns the intercept and derivative
smoothing matrices, in original observation order; this permits reuse of
KernelSmoothing.jl. It must use degree-one Gaussian local regression.

`status=:maxiter` denotes a finite iterate that did not meet the convergence
tolerance, not convergence. No inference guarantee, sign identification, or
support recovery guarantee is attached to this exploratory estimator.
"""
function refine_probe_index(Y::AbstractMatrix,Z::AbstractMatrix,initial::AbstractVector;
        updates::Integer=50,bandwidth_factor::Real=1.0,trim::Real=0.05,
        support=collect(1:size(Z,2)),tolerance::Real=1e-4,
        standardize::Bool=true,smoother=nothing,penalty_method::Symbol=:none,
        penalty_level::Real=0.1)
    n,d=size(Z)
    n>d+2 && size(Y,1)==n && size(Y,2)>0 && length(initial)==d ||
        throw(ArgumentError("invalid MAVE dimensions"))
    all(isfinite,Y) && all(isfinite,Z) && all(isfinite,initial) ||
        throw(ArgumentError("nonfinite MAVE input"))
    updates>=1 && bandwidth_factor>0 && 0<=trim<0.25 && tolerance>0 ||
        throw(ArgumentError("invalid MAVE controls"))
    penalty_method in (:none,:lasso,:adaptive_lasso) && 0<penalty_level<1 ||
        throw(ArgumentError("invalid iterative MAVE penalty"))
    S=sort(unique(collect(support)))
    !isempty(S) && all(j->1<=j<=d,S) || throw(ArgumentError("invalid MAVE support"))
    mu,scale=standardize ? _probe_scale(Y) : (zeros(1,size(Y,2)),ones(1,size(Y,2)))
    V=Matrix((Y .-mu)./scale);b=zeros(d);b[S]=initial[S]
    norm(b)>sqrt(eps(Float64)) || throw(ArgumentError("zero MAVE initialization"))
    b./=norm(b);history=NamedTuple[];status=:maxiter
    penalty_weights=penalty_method==:adaptive_lasso ? 1 ./max.(abs.(b),0.05) : ones(d)
    for it in 1:updates
        P=_mave_system(V,Z,b;bandwidth_factor,trim,smoother)
        lambda=penalty_method==:none ? 0.0 :
            penalty_level*maximum(abs.(P.r[S])./penalty_weights[S])
        penalties=lambda.*penalty_weights
        step=_mave_quadratic(P.Q,P.r;support=S,
                             penalty=penalty_method==:none ? nothing : penalties)
        raw=step.direction
        norm(raw)>sqrt(eps(Float64)) || error("zero MAVE update")
        candidate=normalize(raw);dot(b,candidate)<0 && (candidate=-candidate)
        change=norm(candidate-b)
        before=dot(b,P.Q*b)/2-dot(P.r,b)+dot(penalties,abs.(b))
        after=dot(raw,P.Q*raw)/2-dot(P.r,raw)+dot(penalties,abs.(raw))
        push!(history,(;iteration=it,change,bandwidth=P.h,before,after,ridge=step.ridge,
                       lambda,kkt=step.kkt,direction=copy(candidate)))
        b=candidate
        if updates==1
            status=:one_update;break
        elseif change<tolerance
            status=:converged;break
        end
    end
    (;direction=b,selected=findall(!iszero,b),status,iterations=length(history),history,
      center=vec(mu),scale=vec(scale),bandwidth_factor,trim)
end


# First-order Stein refit on a support chosen without truth. Gaussian covariates
# and a nonzero vector of average probe derivatives identify the index here.
function _mave_stein_refit(Y,Z,selected,reference,method)
    n,d=size(Z);mu,scale=_probe_scale(Y);V=Matrix((Y.-mu)./scale)
    A=Z[:,selected].-sum(Z[:,selected];dims=1)/n
    root=cholesky(Symmetric(A'A/(n-1))).L
    M=root\(A'V/(n-1))
    B=method==:rrr ? M : root'\M
    v=_probe_leading(B*B')
    method==:rrr && (v=root'\v)
    out=zeros(d);out[selected]=v
    out=_orient_index_direction(out)
    # A disjoint new support can be orthogonal to the old direction. Retain the
    # largest-component convention in that case rather than rejecting the fit.
    dot(out,reference)<0 && (out=-out)
    (;direction=out,selected=findall(!iszero,out),status=:stein_refit,iterations=0,
      history=NamedTuple[],center=vec(mu),scale=vec(scale))
end

"""
    select_mave_variables(Y, Z, dense; method=:lasso, level=0.1, refit_updates=50, ...)

Select variables from an input index direction and refit without a penalty
on that support. `:lasso` uses penalty `level * lambda_max`;
`:adaptive_lasso` additionally weights coordinate j by
`1 / max(abs(dense[j]), 0.05)`. Here `lambda_max=max(abs(r[j])/weight[j])`
for the frozen quadratic index criterion at `dense`. `:threshold` retains
`abs(dense[j]) > level * sqrt(log(d)/n)`. These are experimental tuning rules,
not calibrated selection guarantees. The Lasso considers all coordinates.
An empty selection is reported as an error, with no truth-based fallback.
`iterative=true` applies the Lasso rule throughout MAVE starting at `dense`
(which can instead be the Stein initializer), then refits its final support.
`selection_status` and `selection_iterations` describe the penalized stage;
`status` and `iterations` describe the final unpenalized refit.

`refit_method=:rrr` or `:ols` replaces the final MAVE refit by first-order
Stein reduced-rank regression or multivariate OLS on the selected support.
These hybrid options require Gaussian covariates and a nonzero vector of
average probe derivatives. `refit_updates` still controls the penalized stage
when `iterative=true`. The refit has status `:stein_refit` and zero iterations.
A coordinate or support orthogonal to the initializer is not rejected.
"""
function select_mave_variables(Y::AbstractMatrix,Z::AbstractMatrix,dense::AbstractVector;
        method::Symbol=:lasso,level::Real=0.1,refit_updates::Integer=50,
        bandwidth_factor::Real=1.0,trim::Real=0.05,tolerance::Real=1e-4,
        smoother=nothing,iterative::Bool=false,refit_method::Symbol=:mave)
    n,d=size(Z)
    length(dense)==d && norm(dense)>0 && level>0 || throw(ArgumentError("invalid sparse MAVE controls"))
    refit_method in (:mave,:rrr,:ols) || throw(ArgumentError("unknown selected-support refit"))
    b=normalize(dense);kkt=0.0;lambda=0.0
    selection_status=:one_update;selection_iterations=1
    if iterative
        method in (:lasso,:adaptive_lasso) || throw(ArgumentError("iterative selection requires a Lasso rule"))
        selection=refine_probe_index(Y,Z,b;updates=refit_updates,bandwidth_factor,trim,tolerance,
                                     smoother,penalty_method=method,penalty_level=level)
        start=selection.direction;selection_status=selection.status
        selection_iterations=selection.iterations
        lambda=selection.history[end].lambda;kkt=selection.history[end].kkt
    elseif method==:threshold
        start=copy(b);start[abs.(b).<=level*sqrt(log(d)/n)].=0.0
    elseif method in (:lasso,:adaptive_lasso)
        mu,scale=_probe_scale(Y);V=Matrix((Y .-mu)./scale)
        P=_mave_system(V,Z,b;bandwidth_factor,trim,smoother)
        weights=method==:lasso ? ones(d) : 1 ./max.(abs.(b),0.05)
        lambda=level*maximum(abs.(P.r)./weights)
        selection=_mave_quadratic(P.Q,P.r;penalty=lambda.*weights)
        start=selection.direction;kkt=selection.kkt
    else
        throw(ArgumentError("unknown MAVE selection method"))
    end
    selected=findall(!iszero,start)
    isempty(selected) && error("empty MAVE selection")
    fit=if refit_method==:mave
        refine_probe_index(Y,Z,start;updates=refit_updates,bandwidth_factor,trim,
                           support=selected,tolerance,smoother)
    else
        _mave_stein_refit(Y,Z,selected,b,refit_method)
    end
    (;fit...,selection_method=method,level,lambda,kkt,iterative,selection_status,
      selection_iterations,refit_method)
end
