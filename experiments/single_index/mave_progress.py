"""Snapshot diagnostic only. Final primary summaries require all planned shards."""
import csv,statistics,sys
from collections import defaultdict
from pathlib import Path
out=Path(sys.argv[1]);by=defaultdict(list)
for path in sorted(out.glob('part-*.csv')):
 for r in csv.DictReader(path.open()):
  if r['status']=='timing':continue
  if r['method'] in ['stein_screen10','stein_dense20','one_update','full_mave','iterative_adaptive_lasso_refit','holdout_min'] and float(r['bandwidth_factor']) in [0.,1.]:
   by[r['design'],int(r['p']),int(r['n']),r['method']].append(r)
headers=['design','p','n','method','replicates','mean_loss','maxiter','selection_maxiter','mean_missing','mean_extra']
rows=[]
for key,v in sorted(by.items()):
 rows.append([*key,len(v),statistics.fmean(float(r['loss']) for r in v),sum(r['status']=='maxiter' for r in v),sum(r.get('selection_status')=='maxiter' for r in v),statistics.fmean(float(r['missing']) for r in v),statistics.fmean(float(r['extra']) for r in v)])
with (out/'progress-snapshot.csv').open('w') as f:
 w=csv.writer(f);w.writerow(headers);w.writerows(rows)
print(','.join(headers))
for r in rows:
 if r[1]==40 and r[2] in [400,800]:print(','.join(map(str,r)))
