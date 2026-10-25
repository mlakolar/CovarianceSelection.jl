"""Render the prespecified index comparison and small software canary."""
from pathlib import Path

def atom(v):
    if isinstance(v,bool): return str(v).lower()
    if isinstance(v,str): return '"'+v+'"'
    if isinstance(v,list): return '['+', '.join(atom(x) for x in v)+']'
    return str(v)
def save(path, header, scenarios):
    text=''.join(f'{k} = {atom(v)}\n' for k,v in header.items())
    for c in scenarios:
        text+='\n[[scenario]]\n'+''.join(f'{k} = {atom(v)}\n' for k,v in c.items())
    Path(path).write_text(text)

def index_plans(root):
    root=Path(root)
    specs=[(1200,40,20,True),(2400,40,20,True),(2400,80,50,True),
           (2400,40,20,False),(4800,40,20,True),(4800,80,50,True),(9600,80,50,True)]
    scenarios=[]
    for i,(n,p,d,equal) in enumerate(specs,1):
        scenarios.append(dict(id=i,name=f'n{n}_p{p}_d{d}'+('' if equal else '_unequal'),n=n,p=p,d=d,
            index_sparsity=4,index_equal=equal,design='gaussian',graph='fixed_chain',noise='gaussian',
            index_methods=['mean_soft','mean_hard','mean_ols','mom_soft','mom_hard','mom_ols','oracle_support_ols'],
            index_multipliers=[.75,1.,1.25,1.5,2.],replicates=1000))
    save(root/'gaussian_fixed_index.toml',dict(stage='index',seed_base=811000000,shard_size=20,penalty_probability=.1),scenarios)
    canary=[dict(scenarios[0],replicates=2)]
    save(root/'gaussian_fixed_index_canary.toml',dict(stage='index',seed_base=809000000,shard_size=2,penalty_probability=.1),canary)

if __name__=='__main__': index_plans(Path(__file__).parent)
