"""Validate paired records and produce compact reporting inputs under Slurm."""
import csv,json,math,statistics,sys,tomllib
from pathlib import Path
from collections import defaultdict
planpath,out=Path(sys.argv[1]),Path(sys.argv[2])
plan=tomllib.loads(planpath.read_text());groups=defaultdict(list);seen=set();failures=[];timings=[]
scenarios={s['id']:s for s in plan['scenario']}
def configs(s):
    yield ('oracle_probe','ols',20,0.)
    yield ('oracle_probe','true_support_ols',4,0.)
    for bank in (('projections','mixed','complete') if s['p']<=40 else ('projections','mixed')):
        for method in ('ols','moment','rrr'):
            for k in (10,20):yield bank,method,k,1.
        for rho in (.05,.2,.5):
            for method in ('weighted_ols','weighted_rrr'):
                for k in (10,20):yield bank,method,k,rho
        yield bank,'learned_crossfit',20,.2
        yield bank,'setup',0,0.
expected={(s['id'],rep,*c) for s in plan['scenario'] for rep in range(1,s['replicates']+1) for c in configs(s)}
for path in sorted(out.glob('part-*.csv')):
    for r in csv.DictReader(path.open()):
        sid,rep,k=int(r['scenario']),int(r['replicate']),int(r['k']);rho=float(r['rho'])
        key=(sid,rep,r['bank'],r['method'],k,rho)
        assert key not in seen,('duplicate',key)
        seen.add(key)
        if r['status']=='timing':timings.append(float(r['seconds']));continue
        if r['status']!='ok':failures.append(r)
        for f in ('loss','prediction','delta','bandwidth_ratio','weight_l1_small','weight_l1','weight_l1_large','seconds'):r[f]=float(r[f])
        assert 0<=r['loss']<=1 and 0<=r['prediction']<=2
        if r['status']=='ok':assert abs(r['delta']**2-r['prediction'])<1e-8
        r['replicate']=rep;groups[(sid,r['bank'],r['method'],k,rho)].append(r)
assert seen==expected,('missing',len(expected-seen),'extra',len(seen-expected))
def write(path,rows):
    if rows:
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def q90(v):return sorted(v)[math.ceil(.9*len(v))-1]
summary=[];paired=[]
for key,rows in sorted(groups.items()):
    sid,bank,method,k,rho=key;s=scenarios[sid];loss=[r['loss'] for r in rows]
    row=dict(scenario=sid,n=s['n'],p=s['p'],design=s['design'],signal=s['signal'],beta_kind=s['beta_kind'],bank=bank,method=method,k=k,rho=rho,repetitions=len(rows),failures=sum(r['status']!='ok' for r in rows),mean_loss=statistics.mean(loss),mcse_loss=statistics.stdev(loss)/math.sqrt(len(loss)) if len(loss)>1 else 0.,median_loss=statistics.median(loss),q90_loss=q90(loss))
    for f in ('prediction','delta','bandwidth_ratio','weight_l1_small','weight_l1','weight_l1_large','seconds'):
        values=[r[f] for r in rows if math.isfinite(r[f])];row['mean_'+f]=statistics.mean(values) if values else float('nan')
    summary.append(row)
    if bank=='oracle_probe':continue
    baseline=('ols' if method=='weighted_ols' else 'rrr')
    if method in ('weighted_ols','weighted_rrr','learned_crossfit','rrr','moment'):
        if method in ('rrr','moment'):baseline='ols'
        bkey=(sid,bank,baseline,k,1.)
        br={r['replicate']:r for r in groups[bkey]};diff=[r['loss']-br[r['replicate']]['loss'] for r in rows]
        paired.append(dict(scenario=sid,n=s['n'],p=s['p'],design=s['design'],signal=s['signal'],beta_kind=s['beta_kind'],bank=bank,method=method,k=k,rho=rho,baseline=baseline,mean_difference=statistics.mean(diff),mcse=statistics.stdev(diff)/math.sqrt(len(diff)) if len(diff)>1 else 0.,fraction_better=sum(x<0 for x in diff)/len(diff),repetitions=len(diff)))
write(out/'summary.csv',summary);write(out/'paired.csv',paired);write(out/'failures.csv',failures)
receipt=dict(stage=plan['stage'],datasets=sum(s['replicates'] for s in plan['scenario']),rows=len(seen),expected_rows=len(expected),duplicates=0,failures=len(failures),groups=len(summary),setup_seconds=sum(timings),source='inputs-*.sha256',failed_loss_convention=1)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
