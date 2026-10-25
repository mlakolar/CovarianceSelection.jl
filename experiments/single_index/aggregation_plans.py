"""Generate frozen paired designs; does not perform numerical experiments."""
from pathlib import Path
import json
root=Path(__file__).parent
scenarios=[]
def add(p,n,signal=1.,beta='equal'):
    for design in ('chain','pairs'):
        scenarios.append(dict(id=len(scenarios)+1,p=p,n=n,design=design,signal=signal,beta_kind=beta))
for p in (20,40,160):
    for n in (200,800,3200):add(p,n)
for n in (400,1600):add(40,n)
for signal in (.5,1.5):add(40,800,signal)
add(40,800,beta='mixed')
for stage,reps,seed in [('canary',1,729701),('development',100,8729701),('confirmation',500,18729701)]:
    selected=scenarios if stage!='canary' else [scenarios[i] for i in (0,1,10,11,16,17)]
    lines=[f'seed = {seed}',f'shard_size = {1 if stage=="canary" else 10}',f'stage = "{stage}"']
    for s in selected:
        lines+=['','[[scenario]]']+[f'{k} = {json.dumps(v)}' for k,v in s.items()]+[f'replicates = {reps}']
    (root/f'aggregation_{stage}.toml').write_text('\n'.join(lines)+'\n')
