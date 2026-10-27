# Naming migration

Descriptive names use snake_case. Types and modules use UpperCamelCase.
Mathematical identifiers and established numerical abbreviations remain unchanged.
Old function names are deprecated bindings to the same generic function.

| Previous | Canonical |
|---|---|
| `fusedGraphicalLasso` | `fused_graphical_lasso` |
| `fusedGraphicalLasso!` | `fused_graphical_lasso!` |
| `fusedNeighborhoodSelection` | `fused_neighborhood_selection` |
| `groupGraphicalLasso` | `group_graphical_lasso` |
| `groupGraphicalLasso!` | `group_graphical_lasso!` |
| `kernelCovariances` | `kernel_covariances` |
| `tvCovSel` | `tv_covsel` |
| `differencePrecisionActiveShooting` | `difference_precision_active_shooting` |
| `differencePrecisionActiveShooting!` | `difference_precision_active_shooting!` |
| `differencePrecisionRefit` | `difference_precision_refit` |
| `directDifferenceLoss` | `direct_difference_loss` |
| `diffLoss` | `diff_loss` |
| `differencePrecisionIHT` | `difference_precision_iht` |
| `differenceLatentPrecisionIHT` | `difference_latent_precision_iht` |
| `differenceLatentPrecisionIHT!` | `difference_latent_precision_iht!` |
| `covselLatentIHT!` | `covsel_latent_iht!` |
| `HardThreshold!` | `hard_threshold!` |
| `covSel_objective` | `covsel_objective` |
| `getConfusionMatrix` | `get_confusion_matrix` |
| `generateData` | `generate_data` |
| `generateSparseSPD` | `generate_sparse_spd` |
| `genLatentModel` | `gen_latent_model` |
| `relStrength` | `rel_strength` |
| `printEvery` | `print_every` |
| `rowInd` | `row_ind` |
| `colInd` | `col_ind` |

The explicit-step IHT methods and `IHTOptions` were removed in favour of the
automatic-step methods, which choose steps by backtracking. Use
`ProximalOptions` for stopping criteria and callbacks.

Package versions and UUIDs are unchanged. Develop against sibling master checkouts.
