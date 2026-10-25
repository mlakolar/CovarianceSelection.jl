"""Generate prespecified coverage plans. No numerical fitting here."""
import itertools,json
from pathlib import Path
root=Path(__file__).parent
families={'opposed_scale':(.6,'mixed'),'noncommuting':(.6,'complete'),'correlation':(.55,'complete')}

def scenarios(reps,ns=(200,400,800,1600),ratio=1):
    return [dict(id=i+1,design=f,amplitude=a,bank=b,p=40,d=20,n=n,ntrain=n*ratio,replicates=reps) for i,(f,(a,b),n) in enumerate((f,ab,n) for f,ab in families.items() for n in ns)]

def write(stage,ss,configs,indices,seed,shard=10):
    lines=[f'stage = {json.dumps(stage)}',f'seed = {seed}',f'shard_size = {shard}',f'points = [-0.5, 0.0, 0.5]',f'indices = {json.dumps(indices)}']
    for tag,rows in [('scenario',ss),('candidate',configs)]:
        for row in rows:lines+=['',f'[[{tag}]]']+[f'{k} = {json.dumps(v)}' for k,v in row.items()]
    (root/f'promising_coverage_{stage}.toml').write_text('\n'.join(lines)+'\n')

configs=[dict(id=i+1,C=C,hfactor=h,pilot_ratio=r,adjustment='basic') for i,(C,h,r) in enumerate(itertools.product((.75,1.1,1.25),(1.,1.5,2.),(1.,2.,3.)))]
if __name__=='__main__':
    write('canary',scenarios(2,(200,800)),configs,['known','estimated'],81736511,1)
    write('grid',scenarios(300),configs,['known'],91736511,10)
