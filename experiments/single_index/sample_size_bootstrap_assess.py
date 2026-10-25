"""Summarize paired bootstrap intervals without dropping failed replications."""
import csv,json,math,statistics,sys,tomllib
from collections import defaultdict
from pathlib import Path
from gaussian_fixed_assess import finite_mean,wilson

def assess(planpath,out):
    plan=tomllib.loads(Path(planpath).read_text()); out=Path(out)
    scenarios={c['name']:c for c in plan['scenario']}
    shards=sorted(out.glob('shard-*.csv'))
    assert len(shards)==sum(math.ceil(c['replicates']/plan['shard_size']) for c in scenarios.values())
    groups=defaultdict(list); seen=set(); rows=0
    for p in shards:
        with p.open() as f:
            for r in csv.DictReader(f):
                key=tuple(r[k] for k in ('scenario','index_method','interval_method','entry_a','entry_b','point'))
                identity=key+(r['rep'],)
                assert identity not in seen,identity
                seen.add(identity);groups[key].append(r);rows+=1
    results=[]
    for key,rr in sorted(groups.items()):
        n=len(rr); assert n==scenarios[key[0]]['replicates']
        k=sum(int(r['covered']) for r in rr);lo,hi=wilson(k,n)
        valid=[r for r in rr if not r['failure']]
        result=dict(zip(('scenario','index_method','interval_method','entry_a','entry_b','point'),key))
        result.update(n=n,failures=n-len(valid),coverage=k/n,coverage_lower=lo,coverage_upper=hi,
            coverage_among_available=k/len(valid) if valid else float('nan'),
            width=finite_mean(r['width'] for r in rr),
            mean_bootstrap_failures=finite_mean(r['bootstrap_failures'] for r in rr),
            bias=finite_mean(float(r['estimate'])-float(r['truth']) for r in rr),
            empirical_sd=statistics.stdev(float(r['estimate']) for r in rr if math.isfinite(float(r['estimate'])))
                if sum(math.isfinite(float(r['estimate'])) for r in rr)>1 else float('nan'),
            rms_se=math.sqrt(finite_mean(float(r['se'])**2 for r in rr)),
            index_error2=finite_mean(r['index_error2'] for r in rr))
        results.append(result)
    expected=sum(c['replicates']*len(c['index_methods'])*3*len(c['entries'])*c['grid_n'] for c in scenarios.values())
    assert rows==expected,(rows,expected)
    with (out/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),lineterminator='\n');w.writeheader();w.writerows(results)
    receipt=dict(rows=rows,expected_rows=expected,shards=len(shards),groups=len(groups),
        duplicate_identities=0,failures_retained=True,bootstrap_draws=plan['bootstrap_draws'],
        invalid_draw_policy='Any invalid pivot makes that entry/location bootstrap interval unavailable; count as noncoverage.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': assess(sys.argv[1],sys.argv[2])
