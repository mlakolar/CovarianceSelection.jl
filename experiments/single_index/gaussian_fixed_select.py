"""Freeze a feasible index method using development direction error and a gate.
Run under Slurm after gaussian_fixed_assess.py. No inference coverage is read.
"""
import csv, json, math, sys, tomllib
from collections import defaultdict
from pathlib import Path

def freeze(indexplan, resultdir, outpath):
    plan=tomllib.loads(Path(indexplan).read_text());root=Path(resultdir)
    rows=list(csv.DictReader((root/'summary.csv').open()))
    candidates=defaultdict(list)
    byname={c['name']:c for c in plan['scenario']}
    for r in rows:
        if r['method'].startswith('oracle'): continue
        if byname[r['scenario']]['index_equal']:
            candidates[(r['method'],float(r['index_multiplier']))].append(r)
    ranking=[]
    for (method,m),rr in candidates.items():
        passed=[r for r in rr if float(r['exact'])>=.99 and float(r['contains'])>=.995 and int(r['failures'])==0]
        if passed:
            ranking.append(dict(method=method,index_multiplier=m,passed_scenarios=[r['scenario'] for r in passed],
                mean_error2=sum(float(r['error2']) for r in passed)/len(passed)))
    if not ranking:
        decision=dict(status='gate_failed',reason='No feasible candidate meets support gate; no second-stage submission permitted.')
    else:
        # Prefer one estimator that can advance the most prespecified scenarios;
        # compare its direction risk only within its passing scenarios.
        ranking.sort(key=lambda x:(-len(x['passed_scenarios']),x['mean_error2'],x['method'],x['index_multiplier']))
        best=ranking[0]
        passed=[byname[name] for name in best['passed_scenarios']]
        # At most two equal-direction regimes: smallest and largest passing n.
        chosen=[min(passed,key=lambda c:(c['n'],c['p'],c['d']))]
        large=max(passed,key=lambda c:(c['n'],c['p'],c['d']))
        if large['name']!=chosen[0]['name']: chosen.append(large)
        decision=dict(status='gate_passed',selected=best,coverage_scenarios=chosen,
                      criterion='Maximum passing equal-direction scenarios, then minimum mean direction squared error across passing scenarios.',
                      gate=dict(exact=.99,contains=.995,failures=0),ranking=ranking)
    Path(outpath).write_text(json.dumps(decision,indent=2)+'\n')
    print(json.dumps(decision,indent=2))

if __name__=='__main__': freeze(*sys.argv[1:])
