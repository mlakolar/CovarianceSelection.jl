"""Render the next band stage after independently assessing pointwise results."""
import json,sys
from pathlib import Path
from gaussian_fixed_plans import save
freeze=json.loads(Path(sys.argv[1]).read_text());out=Path(sys.argv[2])
assert freeze['status']=='gate_passed'
sel=freeze['selected'];ss=[]
for i,c in enumerate(freeze['coverage_scenarios'],1):
    q=dict(c,id=i,replicates=1000,index_methods=['known',sel['method']],
        index_multipliers=[sel['index_multiplier']],constants=[1.1],h_factors=[1.0],pilot_factors=[2.0],
        variants=['union_refit','corrected'],grid_min=-.75,grid_max=.75,grid_n=21,draws=1999)
    ss.append(q)
# True-nuisance controls share the exact datasets, but do not use an estimated
# index: they are genuine oracle controls for the final smoother and bootstrap.
for q in list(ss):ss.append(dict(q,index_methods=['known'],variants=['oracle_nuisance']))
save(out/'gaussian_fixed_bands.toml',dict(stage='bands',seed_base=877000000,shard_size=10,penalty_probability=.1),ss)
q=dict(ss[0],replicates=200,index_methods=[sel['method']],variants=['union_refit'],grid_n=41)
save(out/'gaussian_fixed_grid_refinement.toml',dict(stage='bands',seed_base=877000000,shard_size=10,penalty_probability=.1),[q])
