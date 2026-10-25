"""Stream shards; retain failed fits and report unconditional coverage."""
import csv,json,math,sys,tomllib
from collections import defaultdict
from pathlib import Path

def report(planpath,out):
    plan=tomllib.loads(Path(planpath).read_text());out=Path(out)
    keys=['scenario','design','n','ntrain','index','config','C','hfactor','pilot_ratio','adjustment','point','target']
    sums=defaultdict(lambda:defaultdict(float));seen=set();failures=[];rows=0
    for path in sorted(out.glob('part-*.csv')):
        for r in csv.DictReader(path.open()):
            ident=tuple(r[k] for k in ('scenario','replicate','index','config','point','target'))
            if ident in seen:raise ValueError(f'duplicate {ident}')
            seen.add(ident);rows+=1
            key=tuple(r[k] for k in keys);s=sums[key];s['replicates']+=1
            ok=r['status']=='ok';s['successful']+=ok;s['covered']+=int(r['covered'])
            s['index_loss']+=float(r['index_loss']);s['orientation_flip']+=int(r['orientation_flip'])
            s['missing']+=int(r['missing']);s['extra']+=int(r['extra'])
            if not ok:
                failures.append({k:r[k] for k in ('scenario','replicate','index','config','point','target','reason')});continue
            estimate,se,truth=map(float,(r['estimate'],r['se'],r['truth']))
            if int(r['covered'])!=int(abs(estimate-truth)<=1.959963984540054*se):raise ValueError('coverage mismatch')
            err=estimate-truth;z=err/se
            for k,v in [('error',err),('error2',err*err),('se',se),('se2',se*se),('z',z),('z2',z*z),('effective_n',float(r['effective_n'])),('selected_mean',float(r['selected_mean'])),('penalty',float(r['penalty']))]:s[k]+=v
    expected=sum(s['replicates']*len(plan['indices'])*sum(c.get('design',s['design'])==s['design'] for c in plan['candidate'])*len(plan['points'])*4 for s in plan['scenario'])
    if rows!=expected:raise ValueError(f'incomplete {rows}/{expected}')
    summary=[]
    for key,s in sorted(sums.items()):
        r=dict(zip(keys,key));n=int(s['replicates']);m=int(s['successful'])
        r.update(replicates=n,successful=m,availability=m/n,coverage=s['covered']/n,conditional_coverage=s['covered']/m if m else math.nan)
        r.update({k:s[k]/n for k in ('index_loss','orientation_flip','missing','extra')})
        r.update({k:s[k]/m if m else math.nan for k in ('se','effective_n','selected_mean','penalty')})
        r['bias']=s['error']/m if m else math.nan
        r['empirical_sd']=math.sqrt(max(0,(s['error2']-s['error']**2/m)/(m-1))) if m>1 else math.nan
        r['rms_se']=math.sqrt(s['se2']/m) if m else math.nan
        r['length']=2*1.959963984540054*r['se']
        r['z_mean']=s['z']/m if m else math.nan
        r['z_sd']=math.sqrt(max(0,(s['z2']-s['z']**2/m)/(m-1))) if m>1 else math.nan
        summary.append(r)
    with (out/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    with (out/'failures.csv').open('w') as f:
        fields=['scenario','replicate','index','config','point','target','reason'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(failures)
    receipt=dict(stage=plan['stage'],rows=rows,expected=expected,unique=len(seen),failures=len(failures),summary_rows=len(summary),coverage_recomputed=True,failed_fits_count_as_noncoverage=True)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':report(*sys.argv[1:])
