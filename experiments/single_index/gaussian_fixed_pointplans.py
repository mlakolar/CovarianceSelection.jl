"""Render second-stage plans only after a saved index-support gate passes."""
import json,sys
from pathlib import Path
from gaussian_fixed_plans import save

if __name__=='__main__':
    freeze=json.loads(Path(sys.argv[1]).read_text());out=Path(sys.argv[2])
    assert freeze['status']=='gate_passed','Index support gate failed; coverage study is blocked.'
    sel=freeze['selected']; scenarios=[]
    for i,c in enumerate(freeze['coverage_scenarios'],1):
        q=dict(c,id=i,replicates=200,index_methods=['known',sel['method']],index_multipliers=[sel['index_multiplier']],
            constants=[.75,1.,1.1,1.25],h_factors=[.6,1.,1.4],pilot_factors=[1.5,2.,3.],
            variants=['plain','union_refit','corrected'],grid_min=-.5,grid_max=.5,grid_n=3,draws=0)
        scenarios.append(q)
    save(out/'gaussian_fixed_point_development.toml',dict(stage='pointwise',seed_base=829000000,shard_size=5,penalty_probability=.1),scenarios)
    q=dict(scenarios[0],replicates=1,constants=[1.1],h_factors=[1.],pilot_factors=[2.],
           variants=['plain','union_refit','corrected','oracle_nuisance'],draws=199)
    save(out/'gaussian_fixed_point_canary.toml',dict(stage='pointwise',seed_base=827000000,shard_size=1,penalty_probability=.1),[q])
