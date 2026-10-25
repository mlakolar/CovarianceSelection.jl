"""Aggregate every outcome, validate counts, and freeze a development shortlist."""
import csv,json,math,statistics,sys,tomllib
from collections import defaultdict,Counter
from pathlib import Path

def write_csv(path,rows):
    if not rows: return
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader();w.writerows(rows)

def config_id(r):
    return f"{r['bank']}-{r['L']}-{r['method']}-{r['k']}"

def report(planfile,out):
    out=Path(out)
    with open(planfile,'rb') as f: plan=tomllib.load(f)
    scenarios={str(s['id']):s for s in plan['scenario']}
    groups=defaultdict(list); identities=set();counts=Counter(); choices=Counter()
    for path in sorted(out.glob('part-*.csv')):
        with path.open() as f:
            for r in csv.DictReader(f):
                s=scenarios[r['scenario']]
                assert 1<=int(r['replicate'])<=s['replicates']
                assert int(r['n'])==s['n'] and int(r['p'])==s['p'] and r['design']==s['design']
                label=r['selected_id'] if r['selector']=='fixed' else r['selector']
                ident=(r['scenario'],r['replicate'],r['selector'],label)
                if ident in identities: raise ValueError(f'duplicate {ident}')
                identities.add(ident);counts[r['scenario']]+=1
                groups[(r['scenario'],r['selector'],label)].append(r)
                if r['selector']=='validation': choices[(r['scenario'],r['selected_id'])]+=1
    nc=len(plan.get('candidate',[])) or 75
    per_rep=nc+(2 if plan.get('validate',False) else 0)
    expected=sum(s['replicates']*per_rep for s in scenarios.values())
    assert len(identities)==expected,(len(identities),expected)
    for sid,s in scenarios.items(): assert counts[sid]==s['replicates']*per_rep
    summary=[]
    for (sid,selector,label),rr in sorted(groups.items()):
        s=scenarios[sid]; valid=[r for r in rr if r['status']=='ok']
        assert len(rr)==s['replicates']
        losses=sorted(float(r['loss']) if r['status']=='ok' else 1.0 for r in rr)
        assert all(math.isfinite(x) and 0<=x<=1 for x in losses)
        pred=[float(r['prediction']) if r['status']=='ok' else 2.0 for r in rr]
        row=dict(scenario=sid,n=s['n'],p=s['p'],design=s['design'],signal=s.get('signal',1.0),
                 beta_kind=s.get('beta_kind','equal'),even=s.get('even',False),selector=selector,
                 id=label,replicates=len(rr),failures=len(rr)-len(valid),
                 mean_loss=statistics.mean(losses),median_loss=statistics.median(losses),
                 q90_loss=losses[math.ceil(.9*len(losses))-1],
                 se_loss=statistics.stdev(losses)/math.sqrt(len(losses)) if len(losses)>1 else 0,
                 prediction=statistics.mean(pred),prob_loss_below_01=sum(x<.1 for x in losses)/len(losses),
                 prob_loss_below_025=sum(x<.25 for x in losses)/len(losses),
                 mean_missing=statistics.mean(float(r['missing']) for r in valid) if valid else math.nan,
                 mean_extra=statistics.mean(float(r['extra']) for r in valid) if valid else math.nan,
                 mean_seconds=statistics.mean(float(r['seconds']) for r in rr))
        summary.append(row)
    write_csv(out/'summary.csv',summary)
    write_csv(out/'selection_counts.csv',[dict(scenario=s,selected_id=c,count=n) for (s,c),n in sorted(choices.items())])
    receipt=dict(rows=len(identities),expected_rows=expected,duplicate_identities=0,
                 scenarios=len(scenarios),configs=nc,failures=sum(r['failures'] for r in summary),
                 all_outcomes_retained=True,plan=str(planfile))
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if not plan.get('candidate'):
        # Equal weight across the twelve main scenarios. Even-link diagnostic
        # results cannot select a winner for the first-order-identifiable study.
        byid=defaultdict(list)
        for r in summary:
            if not r['even']: byid[r['id']].append(r['mean_loss'])
        ranked=sorted((statistics.mean(v),k) for k,v in byid.items())
        write_csv(out/'development_ranking.csv',[dict(id=k,mean_loss=x) for x,k in ranked])
        pools=[lambda k:'-ols-' in k,lambda k:'-screen-' in k,
               lambda k:any('-'+m+'-' in k for m in ('moment','split','second','second_split'))]
        chosen=[]
        for pool in pools:
            chosen.extend([k for x,k in ranked if pool(k)][:2])
        chosen=sorted(set(chosen),key=lambda k:next(x for x,c in ranked if c==k))
        cfg=[]
        for k in chosen:
            bank,L,method,support=k.split('-')
            cfg.append(dict(id=k,bank=bank,L=int(L),method=method,k=int(support)))
        (out/'shortlist.json').write_text(json.dumps(dict(
            rule='Two lowest mean losses in each of OLS, screening/refit, and spectral families; '
                 'equal main-scenario weight; even-link cases excluded; no confirmation data used.',
            candidates=cfg,fixed_default=chosen[0]),indent=2)+'\n')
        print('Development shortlist:',chosen)
    print(json.dumps(receipt))

if __name__=='__main__': report(sys.argv[1],sys.argv[2])
