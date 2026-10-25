"""Full paired hybrid report under Slurm; no removal by support or convergence."""
import csv,json,math,statistics,tomllib
from pathlib import Path
from collections import defaultdict
plan=tomllib.loads(Path('experiments/single_index/mave_confirmation.toml').read_text())
paths=sorted(Path('hybrid').glob('part-*.csv'));assert len(paths)==4800,len(paths)
rows=[r for p in paths for r in csv.DictReader(p.open())];assert len(rows)==9600*8,len(rows)
base={(r['scenario'],r['replicate']):float(r['loss']) for r in rows if r['method']=='stein_screen10'}
by=defaultdict(list);seen=set()
for r in rows:
 key=r['scenario'],r['replicate'],r['method'],r['level'];assert key not in seen;seen.add(key)
 loss=float(r['loss']);assert math.isfinite(loss) and -1e-9<=loss<=1+1e-9
 by[r['scenario'],r['method'],r['level']].append(r)
def quantile(x,q):
 y=sorted(x);v=(len(y)-1)*q;l=math.floor(v)
 return y[l]+(v-l)*(y[min(l+1,len(y)-1)]-y[l])
headers=['scenario' ,'design','p','n','method','level','replicates','mean_loss','mcse_loss','paired_delta','paired_mcse','all_true_rate','exact_support_rate','mean_missing','mean_extra','failed','selection_maxiter','mean_seconds','q90_loss','q95_loss','catastrophic_rate']
with Path('hybrid/summary.csv').open('w') as f:
 w=csv.writer(f);w.writerow(headers)
 for (sid,method,level),v in sorted(by.items()):
  assert len(v)==300
  valid=[r for r in v if r['status']!='failed']
  losses=[float(r['loss']) for r in v];d=[float(r['loss'])-base[r['scenario'],r['replicate']] for r in v]
  w.writerow([sid,v[0]['design'],v[0]['p'],v[0]['n'],method,level,len(v),statistics.fmean(losses),statistics.stdev(losses)/math.sqrt(len(v)),statistics.fmean(d),statistics.stdev(d)/math.sqrt(len(v)),statistics.fmean(int(r['all_true']) for r in v),statistics.fmean(int(r['exact_support']) for r in v),statistics.fmean(float(r['missing']) for r in valid) if valid else math.nan,statistics.fmean(float(r['extra']) for r in valid) if valid else math.nan,sum(r['status']=='failed' for r in v),sum(r['selection_status']=='maxiter' for r in v),statistics.fmean(float(r['seconds']) for r in v),quantile(losses,.9),quantile(losses,.95),statistics.fmean(x>.25 for x in losses)])
receipt=dict(datasets=9600,fit_records=len(rows),shards=len(paths),summary_rows=len(by),failed=sum(r['status']=='failed' for r in rows),oracle='oracle_support_rrr is a labeled diagnostic; no feasible fit receives true support.',validation='Full scenario/replicate/method counts, unique keys, finite bounded losses; pairs share primary seeds; failed fits retain loss one.')
Path('hybrid/receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
