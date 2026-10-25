"""Known-index sample ladder; no index-support gate or tuning by coverage."""
from pathlib import Path
from gaussian_fixed_plans import save
root=Path(__file__).parent
ss=[]
for p in (40,160,400):
    for i,n in enumerate((100,200,400,800,1600,3200),1):
        ss.append(dict(id=i,name=f'chain_n{n}_p{p}',n=n,p=p,d=20,index_sparsity=4,index_equal=True,
            design='gaussian',graph='fixed_chain',noise='gaussian',index_methods=['known'],index_multipliers=[1.25],
            replicates=1000,constants=[1.1],h_factors=[1.0],pilot_factors=[2.0],variants=['union_refit','corrected','oracle_nuisance'],
            entries=[[1,1],[1,2],[1,4]],grid_min=-.5,grid_max=.5,grid_n=3,draws=0))
for i,n in enumerate((100,200,400,800,1600,3200),1):
    ss.append(dict(ss[i-1],id=i,name=f'hub6_n{n}_p160',n=n,p=160,graph='fixed_hub',graph_degree=6,
                   entries=[[1,1],[1,2],[2,3]]))
save(root/'sample_size_known.toml',dict(stage='pointwise',seed_base=911000000,shard_size=50,penalty_probability=.1),ss)
canary=[dict(ss[0],replicates=1),dict(ss[-3],replicates=1),dict(ss[17],replicates=1)]
save(root/'sample_size_canary.toml',dict(stage='pointwise',seed_base=907000000,shard_size=1,penalty_probability=.1),canary)
# Estimated-index comparison keeps all samples, failures included. True-support
# OLS is a diagnostic, not a feasible competitor. The original feasible rule
# remains frozen, even when support selection is unreliable.
ss=[]
for i,n in enumerate((200,400,800,1600,3200),1):
    ss.append(dict(id=i,name=f'index_n{n}_p40',n=n,p=40,d=20,index_sparsity=4,index_equal=True,
        design='gaussian',graph='fixed_chain',noise='gaussian',index_methods=['known','mean_ols','oracle_support_ols'],
        index_multipliers=[1.25],replicates=1000,constants=[1.1],h_factors=[1.0],pilot_factors=[2.0],
        variants=['union_refit'],entries=[[1,1],[1,2],[1,4]],grid_min=-.5,grid_max=.5,grid_n=3,draws=0))
save(root/'sample_size_estimated.toml',dict(stage='pointwise',seed_base=929000000,shard_size=50,penalty_probability=.1),ss)
