include("four_design_common.jl")
using Statistics
import KernelSmoothing
const MV=CovarianceSelection.SingleIndexGraph

# Frozen committed KernelSmoothing version uses smoothermatrix; accept its new name too.
function mv_smoother(t,anchors,h)
    f=isdefined(KernelSmoothing,:smoother_matrix) ? KernelSmoothing.smoother_matrix : KernelSmoothing.smoothermatrix
    k=KernelSmoothing.GaussianKernel()
    (f(t,t[anchors];bandwidth=h,degree=1,derivative=0,kernel=k),
     f(t,t[anchors];bandwidth=h,degree=1,derivative=1,kernel=k))
end

function mv_initial(Y,Z,k=10)
    P=ag_prepare(Y,Z)
    # Exact identity-response-weight RRR, without allocating a large identity matrix.
    if k==size(Z,2)
        return normalize(P.root'\ag_direction(P.N*P.N'))
    end
    selected=sort(sortperm(vec(sum(abs2,P.B;dims=2));rev=true)[1:k])
    A=P.Zc[:,selected];root=cholesky(Symmetric(A'A/(size(A,1)-1))).L
    M=root\(A'P.V/(size(A,1)-1))
    b=zeros(size(Z,2));b[selected]=normalize(root'\ag_direction(M*M'))
    b
end

function mv_predict(Y,Z,b,Znew;bandwidth_factor=1.)
    mu,scale=MV._probe_scale(Y)
    t=Z*b;h=bandwidth_factor*size(Z,1)^(-0.2)*std(t)
    f=isdefined(KernelSmoothing,:smoother_matrix) ? KernelSmoothing.smoother_matrix : KernelSmoothing.smoothermatrix
    A=f(t,Znew*b;bandwidth=h,degree=1,derivative=0,kernel=KernelSmoothing.GaussianKernel())
    prediction=A*((Y.-mu)./scale)
    (;prediction,center=mu,scale)
end

# A nonfinite prediction invalidates that candidate, not every other support.
function mv_support_choice(scores,standard_errors)
    eligible=isfinite.(scores).&isfinite.(standard_errors).&(standard_errors.>=0)
    adjusted=ifelse.(eligible,scores,Inf)
    best=argmin(adjusted)
    isfinite(adjusted[best]) || throw(DomainError(scores,"no finite validation candidate"))
    sparse=findfirst(eligible .& (scores .<=scores[best]+standard_errors[best]))
    sparse===nothing && error("finite validation minimum has no one-SE candidate")
    best,sparse,eligible
end

# Honest support selection. Split fixed independently of responses. Selection sees
# only 75% training data. Refit the chosen support on all observations afterwards.
function mv_holdout(Y,Z;seed,bandwidth_factor=1.,updates=40)
    n,d=size(Z);perm=randperm(MersenneTwister(seed),n);nt=floor(Int,0.75n)
    train=perm[1:nt];val=perm[nt+1:end]
    Yt=Y[train,:];Zt=Z[train,:]
    initial=mv_initial(Yt,Zt)
    dense=refine_probe_index(Yt,Zt,initial;updates,bandwidth_factor,smoother=mv_smoother)
    ranking=sortperm(abs.(dense.direction);rev=true)
    sizes=[2,4,6,10,d]
    fits=NamedTuple[];scores=Float64[];standard_errors=Float64[];nonfinite_losses=Int[]
    for k in sizes
        S=sort(ranking[1:k]);start=zeros(d);start[S]=dense.direction[S]
        fit=k==d ? dense : refine_probe_index(Yt,Zt,start;updates,bandwidth_factor,
                                              support=S,smoother=mv_smoother)
        pred=mv_predict(Yt,Zt,fit.direction,Z[val,:];bandwidth_factor)
        losses=vec(mean(abs2.((Y[val,:].-pred.center)./pred.scale.-pred.prediction);dims=2))
        push!(fits,fit);push!(scores,mean(losses));push!(standard_errors,std(losses)/sqrt(length(val)))
        push!(nonfinite_losses,count(!isfinite,losses))
    end
    # Minimum validation error is primary. One-SE is a separate sparse rule.
    best,sparse,eligible=mv_support_choice(scores,standard_errors)
    results=Dict{String,Any}()
    for (name,j) in (("holdout_min",best),("holdout_1se",sparse))
        fit=fits[j]
        full=refine_probe_index(Y,Z,fit.direction;updates,bandwidth_factor,
                                support=fit.selected,smoother=mv_smoother)
        results[name]=(fit=full,k=sizes[j],validation=scores[j],training_status=String(fit.status))
    end
    (;results,scores,sizes,standard_errors,nonfinite_losses,eligible,training_status=String(dense.status))
end
