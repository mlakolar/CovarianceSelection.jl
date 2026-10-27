"""
    CovarianceSelection

Methods for learning high-dimensional Gaussian graphical models:

* CLIME and ACLIME, sparse precision matrices by constrained ℓ1 minimization
  ([`clime`](@ref), [`aclime`](@ref)), with cross validation
  ([`climecv`](@ref), [`aclimecv`](@ref)) and a choice of solvers
  ([`CLIMESolver`](@ref));
* the precision matrix by the scaled lasso of each column, which needs no
  tuning ([`slassoinv`](@ref));
* the graphical lasso by ADMM ([`covsel!`](@ref), [`covselpath`](@ref));
* sparse minus low rank precision matrices of graphical models with latent
  variables by iterative hard thresholding ([`covsel_latent_iht!`](@ref));
* the fused graphical lasso for K groups ([`fused_graphical_lasso`](@ref)) and
  fused neighborhood selection ([`fused_neighborhood_selection`](@ref));
* K precision matrices with a shared sparsity pattern, the group graphical
  lasso ([`group_graphical_lasso`](@ref)), and graphs that vary over time by
  kernel smoothing ([`tv_covsel`](@ref), [`kernel_covariances`](@ref));
* direct estimation of the difference of two precision matrices by coordinate
  descent ([`difference_precision_active_shooting`](@ref),
  [`difference_precision_refit`](@ref)) or iterative hard thresholding, with an
  optional low rank component ([`difference_precision_iht`](@ref),
  [`difference_latent_precision_iht`](@ref)).

The losses and penalties are those of ProximalBase, and the graphical lasso and
the direct difference estimator are solved by the ADMM and the coordinate
descent of ProximalSolvers. The IHT estimators use ProximalSolvers' proximal
gradient methods. Before its first release the package was called CovSel.
"""
module CovarianceSelection

using ProximalBase
using ProximalSolvers: ProximalSolvers, ADMMSolver, CDSolver, ProximalOptions, refit!, solve!
using SparseArrays, LinearAlgebra
using Random: AbstractRNG, default_rng, randperm
using SpecialFunctions: erfcinv

export
  ADMMOptions, ADMMInfo,
  covsel!, covselpath,
  covsel_latent_iht!,

  # difference estimation
  #
  fused_graphical_lasso, fused_graphical_lasso!,   # Danaher et al.
  fused_neighborhood_selection,

  # K precision matrices: shared support, time-varying graphs
  group_graphical_lasso, group_graphical_lasso!, kernel_covariances, tv_covsel,

  # Direct Difference Estimation and Inference
  difference_precision_active_shooting, difference_precision_active_shooting!, difference_precision_refit,
  direct_difference_loss, diff_loss,

  difference_precision_iht,
  difference_latent_precision_iht, difference_latent_precision_iht!,

  # CLIME and ACLIME
  clime, aclime, CLIMEPath, climecv, aclimecv, CLIMECV, bestestimate,
  CLIMESolver, DualSimplex, ParametricSimplex, LPSolver, PDHG,

  # scaled lasso of Sun and Zhang
  slassoinv

include("utils.jl")

include("admmCovarianceSelection.jl")

include("directDifference.jl")
include("automaticIHT.jl")
include("fusedDifference.jl")
include("group_graphical_lasso.jl")

include("clime_simplex.jl")
include("clime.jl")
include("clime_pdhg.jl")
include("clime_cv.jl")
include("slassoinv.jl")

include("generate.jl")
include("evaluation.jl")

include("deprecated_naming.jl")
include("precompile.jl")

end
