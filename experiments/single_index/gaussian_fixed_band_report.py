"""Band summaries and paired 21/41-grid differences from saved results."""
import csv,json,math,statistics,sys
from pathlib import Path
root=Path(sys.argv[1]);out=root/'results'
rows=list(csv.DictReader((out/'summary.csv').open()))
bands=[r for r in rows if r['kind'] in ('curve','family','family_component')]
fields=['scenario','index_method','variant','kind','entry_a','entry_b','n','failures','coverage','coverage_lower','coverage_upper','width','index_exact']
with (out/'bands.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(bands)
print(json.dumps([{k:r[k] for k in fields} for r in bands if r['kind']=='family'],indent=2))

def read_paired(directory):
    selected={}
    for path in sorted(directory.glob('shard-*.csv')):
        with path.open() as f:
            for r in csv.DictReader(f):
                if r['scenario']!='n4800_p40_d20' or r['index_method']!='mean_ols' or r['variant']!='union_refit' or int(r['rep'])>200 or r['kind'] not in ('curve','family'): continue
                key=(r['kind'],r['entry_a'],r['entry_b'],int(r['rep']))
                assert key not in selected
                selected[key]=r
    return selected
if (root/'grid41/receipt.json').exists():
    base=read_paired(out);ref=read_paired(root/'grid41');assert base.keys()==ref.keys() and len(base)==800
    comparisons=[]
    for kind,a,b in sorted(set(k[:3] for k in base)):
        keys=[k for k in base if k[:3]==(kind,a,b)]
        diffs=[int(ref[k]['covered'])-int(base[k]['covered']) for k in keys]
        comparisons.append(dict(kind=kind,entry_a=int(a),entry_b=int(b),n=len(keys),
            coverage21=statistics.mean(int(base[k]['covered']) for k in keys),
            coverage41=statistics.mean(int(ref[k]['covered']) for k in keys),
            difference=statistics.mean(diffs),paired_mcse=statistics.stdev(diffs)/math.sqrt(len(diffs)),
            width21=statistics.mean(float(base[k]['width']) for k in keys),
            width41=statistics.mean(float(ref[k]['width']) for k in keys)))
    (out/'paired_grid_refinement.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    print(json.dumps(comparisons,indent=2))
