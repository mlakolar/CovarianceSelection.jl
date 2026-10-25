"""Freeze settings and fixed estimator rules; no numerical computation."""
from pathlib import Path
import json
root=Path(__file__).parent
families={'correlation':[.15,.35,.55],'opposed_scale':[.2,.4,.6],'fixed_diagonal':[.2,.4,.6],'noncommuting':[.2,.4,.6]}
scenarios=[]
for design,amps in families.items():
    for a in amps:
        for p in (20,40):
            for d in (4,20):
                for n in (100,200,400,800,1600):
                    scenarios.append(dict(id=len(scenarios)+1,design=design,amplitude=a,p=p,d=d,n=n))
configs=[]
for d in (4,20):
    ks=[d] if d==4 else [10,20]
    for bank in ('global_trace','global_products'):
        for k in ks:configs.append(dict(d=d,bank=bank,method='ols',k=k,rho=1.))
    configs += [dict(d=d,bank='oracle_probe',method='ols',k=d,rho=0.),dict(d=d,bank='oracle_probe',method='true_support_ols',k=4,rho=0.)]
    for bank in ('projections','mixed','complete'):
        for method in ('ols','moment','rrr'):
            for k in ks:configs.append(dict(d=d,bank=bank,method=method,k=k,rho=1.))
        for rho in (.05,.2,.5):
            for method in ('weighted_ols','weighted_rrr'):
                for k in ks:configs.append(dict(d=d,bank=bank,method=method,k=k,rho=rho))
        configs.append(dict(d=d,bank=bank,method='learned_crossfit',k=d,rho=.2))
        configs.append(dict(d=d,bank=bank,method='split',k=d,rho=1.))
        if bank=='mixed':configs.append(dict(d=d,bank=bank,method='second',k=d,rho=1.))
for stage,reps,seed in [('canary',1,39473811),('comparison',300,59470271)]:
    selected=scenarios if stage=='comparison' else [s for s in scenarios if s['amplitude']==families[s['design']][-1] and ((s['p'],s['d'],s['n']) in ((20,4,100),(40,20,1600)))]
    lines=[f'stage = "{stage}"',f'seed = {seed}',f'shard_size = {1 if stage=="canary" else 25}']
    for s in selected:lines+=['','[[scenario]]']+[f'{k} = {json.dumps(v)}' for k,v in s.items()]+[f'replicates = {reps}']
    for c in configs:lines+=['','[[candidate]]']+[f'{k} = {json.dumps(v)}' for k,v in c.items()]
    (root/f'four_design_{stage}.toml').write_text('\n'.join(lines)+'\n')
