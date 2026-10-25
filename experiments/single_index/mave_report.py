"""Aggregate paired index simulations. Run under Slurm, including validation."""
import csv,json,math,statistics,sys,tomllib
from collections import defaultdict
from pathlib import Path
plan=tomllib.loads(Path(sys.argv[1]).read_text());out=Path(sys.argv[2])
scenarios={s['id']:s for s in plan['scenario']}
expected_shards=sum(math.ceil(s['replicates']/plan['shard_size']) for s in plan['scenario'])
paths=sorted(out.glob('part-*.csv'))
assert len(paths)==expected_shards,(len(paths),expected_shards)
rows=[]
for path in paths:
    with path.open() as f:rows.extend(csv.DictReader(f))
fits=[r for r in rows if r['status']!='timing']
seen=set();by=defaultdict(list);baseline={}
for r in fits:
    sid,rep=int(r['scenario']),int(r['replicate']);s=scenarios[sid]
    assert 1<=rep<=s['replicates']
    key=(sid,rep,r['method'],r['bandwidth_factor'],r['level'])
    assert key not in seen,key
    seen.add(key)
    loss=float(r['loss']);assert math.isfinite(loss) and -1e-9<=loss<=1+1e-9
    if r['method']=='stein_screen10':baseline[sid,rep]=loss
    by[(sid,r['method'],r['bandwidth_factor'],'0.0' if r['method'].startswith('holdout_') else r['level'])].append(r)
expected_methods=2+len(plan['bandwidth_factors'])*(2+len(plan['threshold_levels'])+2*len(plan['lasso_levels']))+2*int(plan['holdout'])+2*int(plan.get('iterative_selection',False))
assert len(fits)==sum(s['replicates'] for s in scenarios.values())*expected_methods,(len(fits),expected_methods)
def mean(x):return statistics.fmean(x)
def se(x):return statistics.stdev(x)/math.sqrt(len(x)) if len(x)>1 else 0.
def quantile(x,q):
    y=sorted(x);v=(len(y)-1)*q;l=math.floor(v);return y[l]+(v-l)*(y[min(l+1,len(y)-1)]-y[l])
headers=['scenario','design','p','d','n','bank','method','bandwidth_factor','level','replicates','mean_loss','mcse_loss','median_loss','q90_loss','q95_loss','catastrophic_rate','paired_delta','paired_mcse','improved_fraction','mean_missing','mean_extra','all_true_rate','exact_support_rate','failures','maxiter','converged','mean_iterations','mean_seconds','mean_final_change','selection_maxiter']
summary=[]
for (sid,method,h,level),group in sorted(by.items()):
    s=scenarios[sid];assert len(group)==s['replicates'],(sid,method,len(group))
    losses=[float(r['loss']) for r in group]
    differences=[float(r['loss'])-baseline[sid,int(r['replicate'])] for r in group]
    valid=[r for r in group if r['status']!='failed']
    changes=[float(r['change']) for r in valid if math.isfinite(float(r['change']))]
    summary.append([sid,s['design'],s['p'],20,s['n'],s['bank'],method,h,level,len(group),mean(losses),se(losses),quantile(losses,.5),quantile(losses,.9),quantile(losses,.95),mean([x>.25 for x in losses]),mean(differences),se(differences),mean([x<0 for x in differences]),mean([float(r['missing']) for r in valid]) if valid else math.nan,mean([float(r['extra']) for r in valid]) if valid else math.nan,mean([int(r['all_true']) for r in group]),mean([int(r['exact_support']) for r in group]),sum(r['status']=='failed' for r in group),sum(r['status']=='maxiter' for r in group),sum(r['status']=='converged' for r in group),mean([int(r['iterations']) for r in group]),mean([float(r['seconds']) for r in group]),mean(changes) if changes else math.nan,sum(r.get('selection_status')=='maxiter' for r in group)])
with (out/'summary.csv').open('w') as f:
    w=csv.writer(f);w.writerow(headers);w.writerows(summary)
receipt=dict(datasets=sum(s['replicates'] for s in scenarios.values()),fit_records=len(fits),timing_records=len(rows)-len(fits),summary_rows=len(summary),shards=len(paths),failed_fits=sum(r['status']=='failed' for r in fits),maxiter_fits=sum(r['status']=='maxiter' for r in fits),validation='Expected shard, dataset, method and replicate counts; unique keys; finite bounded loss; paired baseline matching. Maxiter iterates retained and labeled separately from convergence.',plan=sys.argv[1])
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
