"""Prespecified independent confirmation of the central fixed-bandwidth rule.
Render before outcomes; submit after inspecting the development surface. No truth enters
estimation on a replication. The constants are fixed between replications.
"""
import json,sys
from pathlib import Path
from gaussian_fixed_plans import save
frozen=json.loads(Path(sys.argv[1]).read_text());out=Path(sys.argv[2])
assert frozen['status']=='gate_passed'
sel=frozen['selected'];ss=[]
for i,c in enumerate(frozen['coverage_scenarios'],1):
    ss.append(dict(c,id=i,replicates=1000,index_methods=['known',sel['method']],
        index_multipliers=[sel['index_multiplier']],constants=[.75,1.1,1.25],h_factors=[1.0],
        pilot_factors=[2.0],variants=['plain','union_refit','corrected'],grid_min=-.5,grid_max=.5,grid_n=3,draws=0))
save(out/'gaussian_fixed_point_confirmation.toml',dict(stage='pointwise',seed_base=853000000,shard_size=20,penalty_probability=.1),ss)
