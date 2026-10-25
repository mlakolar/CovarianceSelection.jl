"""Streaming paired validation/reporting. Execute under Discovery Slurm."""
import csv,json,math,statistics,sys,tomllib
from pathlib import Path
from collections import defaultdict
planpath,out=Path(sys.argv[1]),Path(sys.argv[2]);plan=tomllib.loads(planpath.read_text())
scenarios={s['id']:s for s in plan['scenario']};configs={d:[c for c in plan['candidate'] if c['d']==d] for d in (4,20)}
fields=('prediction','delta','bandwidth_ratio','weight_l1_small','weight_l1','weight_l1_large','missing','extra','seconds')
def configkey(c):return (c['bank'],c['method'],int(c['k']),float(c['rho']))
expected={d:{configkey(c) for c in configs[d]}|{(b,'setup',0,0.) for b in ('projections','mixed','complete')} for d in (4,20)}
assert all(len(configs[d])==len(expected[d])-3 for d in configs)
groups={};pairs=defaultdict(list);seen=set();rows_total=0;failures=[];setup_seconds=0.
def consume(sid,rep,records):
    global rows_total,setup_seconds
    s=scenarios[sid];assert (sid,rep) not in seen;seen.add((sid,rep))
    assert set(records)==expected[s['d']],('config mismatch',sid,rep)
    for key,r in records.items():
        rows_total+=1
        if r['status']=='timing':setup_seconds+=float(r['seconds']);continue
        if r['status']!='ok':failures.append(r.copy())
        values={f:float(r[f]) for f in ('loss',)+fields}
        assert 0<=values['loss']<=1 and 0<=values['prediction']<=2
        assert abs(values['delta']**2-values['prediction'])<1e-7
        gkey=(sid,*key)
        g=groups.setdefault(gkey,dict(loss=[],sums={f:0. for f in fields},counts={f:0 for f in fields},failures=0))
        g['loss'].append(values['loss']);g['failures']+=r['status']!='ok'
        for f in fields:
            if math.isfinite(values[f]) and not(f in ('missing','extra') and r['status']!='ok'):
                g['sums'][f]+=values[f];g['counts'][f]+=1
        bank,method,k,rho=key
        if bank in ('global_trace','global_products','oracle_probe'):continue
        baseline=None
        if method in ('weighted_ols','weighted_rrr'):baseline=method.removeprefix('weighted_')
        elif method in ('moment','rrr'):baseline='ols'
        elif method in ('learned_crossfit','split','second'):baseline='rrr'
        if baseline:
            b=records[(bank,baseline,k,1.)]
            pairs[(sid,bank,method,k,rho,baseline)].append(values['loss']-float(b['loss']))
for path in sorted(out.glob('part-*.csv')):
    batch={}
    for r in csv.DictReader(path.open()):
        sid,rep=int(r['scenario']),int(r['replicate']);s=scenarios[sid]
        assert all(str(r[f])==str(s[f]) for f in ('design','p','d','n'))
        assert float(r['amplitude'])==s['amplitude']
        key=(r['bank'],r['method'],int(r['k']),float(r['rho']))
        records=batch.setdefault((sid,rep),{});assert key not in records;records[key]=r
    for (sid,rep),records in batch.items():consume(sid,rep,records)
expected_datasets={(s['id'],r) for s in plan['scenario'] for r in range(1,s['replicates']+1)}
assert seen==expected_datasets,('dataset count',len(seen),len(expected_datasets))
summary=[];paired=[]
for (sid,bank,method,k,rho),g in sorted(groups.items()):
    s=scenarios[sid];v=g['loss'];n=len(v)
    assert n==s['replicates']
    r=dict(scenario=sid,design=s['design'],amplitude=s['amplitude'],p=s['p'],d=s['d'],n=s['n'],bank=bank,method=method,k=k,rho=rho,repetitions=n,failures=g['failures'],mean_loss=statistics.mean(v),mcse_loss=statistics.stdev(v)/math.sqrt(n) if n>1 else 0.,median_loss=statistics.median(v),q90_loss=sorted(v)[math.ceil(.9*n)-1])
    for f in fields:r['mean_'+f]=g['sums'][f]/g['counts'][f] if g['counts'][f] else float('nan')
    summary.append(r)
for (sid,bank,method,k,rho,baseline),v in sorted(pairs.items()):
    s=scenarios[sid];assert len(v)==s['replicates']
    paired.append(dict(scenario=sid,design=s['design'],amplitude=s['amplitude'],p=s['p'],d=s['d'],n=s['n'],bank=bank,method=method,k=k,rho=rho,baseline=baseline,mean_difference=statistics.mean(v),mcse=statistics.stdev(v)/math.sqrt(len(v)) if len(v)>1 else 0.,fraction_better=sum(x<0 for x in v)/len(v),repetitions=len(v)))
def write(path,rows):
    if rows:
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write(out/'summary.csv',summary);write(out/'paired.csv',paired);write(out/'failures.csv',failures)
expected_rows=sum(s['replicates']*len(expected[s['d']]) for s in plan['scenario'])
assert rows_total==expected_rows
receipt=dict(stage=plan['stage'],scenarios=len(scenarios),datasets=len(seen),rows=rows_total,expected_rows=expected_rows,duplicate_identities=0,failures=len(failures),groups=len(groups),setup_seconds=setup_seconds,failed_loss_convention=1,precision_estimation=False,inference=False)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
