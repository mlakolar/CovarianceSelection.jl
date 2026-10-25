"""Compact stage summaries. Run on Discovery, including all failed fits."""
import csv, json, math, statistics, sys, tomllib
from collections import defaultdict
from pathlib import Path

def finite_mean(values):
    a=[float(x) for x in values if math.isfinite(float(x))]
    return statistics.mean(a) if a else float('nan')

def wilson(k,n):
    z=1.959963984540054
    p=k/n; den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    rad=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-rad,mid+rad

def assess(planpath, root):
    plan=tomllib.loads(Path(planpath).read_text()); root=Path(root)
    index=plan['stage']=='index'; groups=defaultdict(list); seen=set(); rows=0
    expected_shards=sum(math.ceil(c['replicates']/plan['shard_size']) for c in plan['scenario'])
    shards=sorted(root.glob('shard-*.csv'))
    assert len(shards)==expected_shards,(len(shards),expected_shards)
    for path in shards:
        with path.open() as f:
            for row in csv.DictReader(f):
                rows+=1
                if index:
                    key=tuple(row[k] for k in ['scenario','method','index_multiplier'])
                    identity=key+(row['rep'],)
                else:
                    key=tuple(row[k] for k in ['scenario','index_method','index_multiplier','C','h_factor','pilot_factor','variant','kind','entry_a','entry_b','point'])
                    identity=key+(row['rep'],)
                assert identity not in seen,identity
                seen.add(identity); groups[key].append(row)
    results=[]
    for key,rr in sorted(groups.items()):
        n=len(rr)
        if index:
            names=['scenario','method','index_multiplier']
            out=dict(zip(names,key)); exact=sum(int(r['exact']) for r in rr)
            lo,hi=wilson(exact,n)
            out.update(n=n,failures=sum(bool(r['failure']) for r in rr),error2=finite_mean(r['error2'] for r in rr),
                exact=exact/n,exact_lower=lo,exact_upper=hi,contains=finite_mean(r['contains'] for r in rr),
                missed=finite_mean(r['missed'] for r in rr),extra=finite_mean(r['extra'] for r in rr),seconds=finite_mean(r['seconds'] for r in rr))
        else:
            names=['scenario','index_method','index_multiplier','C','h_factor','pilot_factor','variant','kind','entry_a','entry_b','point']
            out=dict(zip(names,key)); covered=sum(int(r['covered']) for r in rr)
            lo,hi=wilson(covered,n)
            estimates=[float(r['estimate']) for r in rr if not r['failure'] and math.isfinite(float(r['estimate']))]
            errors=[float(r['estimate'])-float(r['truth']) for r in rr if not r['failure'] and math.isfinite(float(r['estimate']))]
            out.update(n=n,failures=sum(bool(r['failure']) for r in rr),coverage=covered/n,coverage_lower=lo,coverage_upper=hi,
                width=finite_mean(r['width'] for r in rr),bias=finite_mean(errors),
                empirical_sd=statistics.stdev(estimates) if len(estimates)>1 else float('nan'),
                rms_se=math.sqrt(finite_mean(float(r['se'])**2 for r in rr)),
                mse=finite_mean(x*x for x in errors),index_error2=finite_mean(r['index_error2'] for r in rr),
                index_exact=finite_mean(r['index_exact'] for r in rr),index_missed=finite_mean(r['index_missed'] for r in rr),
                index_extra=finite_mean(r['index_extra'] for r in rr),nuisance_missed=finite_mean(r['nuisance_missed'] for r in rr),
                nuisance_contains=finite_mean(float(r['nuisance_missed'])==0 for r in rr if math.isfinite(float(r['nuisance_missed']))),
                nuisance_extra=finite_mean(r['nuisance_extra'] for r in rr),nuisance_size=finite_mean(r['nuisance_size'] for r in rr),
                effective_n=finite_mean(r['effective_n'] for r in rr),lambda_mean=finite_mean(r['lambda'] for r in rr),
                final_effective_n=finite_mean(r.get('final_effective_n','nan') for r in rr),
                final_local_n=finite_mean(r.get('final_local_n','nan') for r in rr))
        results.append(out)
    # Every setting must retain all requested replications, even failures.
    byname={c['name']:c for c in plan['scenario']}
    assert all(r['n']==byname[r['scenario']]['replicates'] for r in results)
    expected=0
    for c in plan['scenario']:
        modes=sum(1 if m in ('known','oracle_support_ols') else len(c['index_multipliers']) for m in c['index_methods'])
        if index: expected+=c['replicates']*modes
        else: expected+=c['replicates']*modes*len(c['constants'])*len(c['h_factors'])*len(c['pilot_factors'])*len(c['variants'])*(3*c['grid_n']+(7 if c['draws']>=100 else 0))
    assert rows==expected,(rows,expected)
    with (root/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),lineterminator="\n");w.writeheader();w.writerows(results)
    receipt=dict(stage=plan['stage'],rows=rows,expected_rows=expected,shards=len(shards),duplicate_identities=0,
                 failures_retained=True,groups=len(results),plan=str(planpath))
    (root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if index:
        ranked=defaultdict(list)
        for r in results:
            if not r['method'].startswith('oracle'):
                ranked[(r['method'],r['index_multiplier'])].append(r)
        ranking=[]
        for (method,m),rr in ranked.items():
            ranking.append(dict(method=method,index_multiplier=float(m),mean_error2=statistics.mean(r['error2'] for r in rr),
                minimum_exact=min(r['exact'] for r in rr),minimum_contains=min(r['contains'] for r in rr),
                maximum_failures=max(r['failures'] for r in rr),gate_pass=all(r['exact']>=.99 and r['contains']>=.995 and r['failures']==0 for r in rr)))
        ranking.sort(key=lambda r:(not r['gate_pass'],r['mean_error2']))
        (root/'ranking.json').write_text(json.dumps(ranking,indent=2)+'\n')
        print(json.dumps(dict(receipt=receipt,best=ranking[:5]),indent=2))
    else: print(json.dumps(receipt,indent=2))

if __name__=='__main__': assess(sys.argv[1],sys.argv[2])
