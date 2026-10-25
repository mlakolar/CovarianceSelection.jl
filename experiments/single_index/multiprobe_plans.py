"""Write fixed, reproducible stage plans. No fitting runs in this script."""
import json,sys
from pathlib import Path

def write(path,seed,scenarios,candidates=None,validate=False):
    lines=[f'seed = {seed}','shard_size = 100',f'validate = {str(validate).lower()}']
    for c in candidates or []:
        lines+=['','[[candidate]]']+[f'{k} = {json.dumps(v)}' for k,v in c.items()]
    for i,s in enumerate(scenarios,1):
        lines+=['','[[scenario]]',f'id = {i}']+[f'{k} = {json.dumps(v)}' for k,v in s.items()]
    Path(path).write_text('\n'.join(lines)+'\n')

def main(stage,out,shortlist=None):
    if stage=='smoke':
        cases=[('diagonal',0,'ols',20),('complete',0,'ols',20),('mixed',30,'screen',6),
               ('mixed',30,'split',20),('mixed',30,'second',20),('mixed',30,'second_split',20)]
        candidates=[dict(id=f'{b}-{L}-{m}-{k}',bank=b,L=L,method=m,k=k) for b,L,m,k in cases]
        write(out,19012791,[dict(n=100,p=40,design='chain',replicates=2),
                             dict(n=200,p=40,design='pairs',replicates=2)],candidates,True)
    elif stage=='development':
        scenarios=[dict(n=n,p=40,design=d,replicates=300)
                   for d in ('chain','pairs') for n in (100,200,400,800,1600,3200)]
        scenarios += [dict(n=n,p=40,design=d,replicates=200,even=True)
                      for d in ('chain','pairs') for n in (400,800,1600)]
        write(out,27012511,scenarios)
    elif stage=='confirmation':
        candidates=json.loads(Path(shortlist).read_text())['candidates']
        scenarios=[dict(n=n,p=40,design=d,replicates=1000)
                   for d in ('chain','pairs') for n in (100,200,400,800,1600,3200)]
        scenarios += [dict(n=n,p=160,design=d,replicates=1000)
                      for d in ('chain','pairs') for n in (400,800,1600)]
        scenarios += [dict(n=n,p=40,design=d,signal=signal,replicates=1000)
                      for d in ('chain','pairs') for n in (400,1600) for signal in (0.5,1.5)]
        scenarios += [dict(n=n,p=40,design=d,beta_kind='mixed',replicates=1000)
                      for d in ('chain','pairs') for n in (400,1600)]
        write(out,57012631,scenarios,candidates,True)
    else: raise ValueError(stage)

if __name__=='__main__': main(*sys.argv[1:])
