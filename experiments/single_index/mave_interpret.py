"""Compact paired comparisons; run inside a Discovery allocation."""
import csv,json,math,statistics,sys
from collections import defaultdict
from pathlib import Path
root=Path(sys.argv[1]);out=root/'interpretation';out.mkdir(exist_ok=True)
lines=['# MAVE index-only comparisons','','All losses are squared sine direction errors. Numerical failures receive loss one. Finite maxiter fits remain in the primary averages. No precision estimation or inference is included.','']
for stage in ['development','confirmation']:
 p=root/stage/'summary.csv'
 if not p.exists():continue
 rows=list(csv.DictReader(p.open()))
 receipt=json.loads((root/stage/'receipt.json').read_text())
 lines += [f'## {stage.title()}','',f"{receipt['datasets']} datasets; {receipt['fit_records']} fits; {receipt['failed_fits']} failures; {receipt['maxiter_fits']} maxiter fits.",'']
 lines += ['| Design | p | n | Stein10 | Stein20 | One update | Full MAVE | Threshold/refit | Adaptive/refit | Iterative adaptive/refit | Holdout min |','| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
 groups=defaultdict(dict)
 for r in rows:
  m=r['method'];h=float(r['bandwidth_factor']);level=float(r['level'])
  if (m.startswith('stein_') or h==1.) and (m not in ['threshold_refit','lasso_refit','adaptive_lasso_refit'] or level==(1. if m=='threshold_refit' else .1)):
   groups[r['design'],int(r['p']),int(r['n'])][m]=r
 methods=['stein_screen10','stein_dense20','one_update','full_mave','threshold_refit','adaptive_lasso_refit','iterative_adaptive_lasso_refit','holdout_min']
 for (design,p,n),g in sorted(groups.items()):
  lines.append('| '+' | '.join([design,str(p),str(n)]+[f"{float(g[m]['mean_loss']):.5f}" if m in g else 'NA' for m in methods])+' |')
 lines += ['','## Support containment and convergence (p=40, n=800; main bandwidth)','', '| Design | Method | Mean loss | Paired difference (MCSE) | All true | Exact | Extras | Maxiter | Penalized/training maxiter |','| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
 for (design,p,n),g in sorted(groups.items()):
  if p!=40 or n!=800:continue
  for m in methods+['iterative_lasso_refit','holdout_1se']:
   if m not in g:continue
   r=g[m];lines.append('| '+' | '.join([design,m,f"{float(r['mean_loss']):.5f}",f"{float(r['paired_delta']):+.5f} ({float(r['paired_mcse']):.5f})",f"{100*float(r['all_true_rate']):.1f}%",f"{100*float(r['exact_support_rate']):.1f}%",f"{float(r['mean_extra']):.2f}",r['maxiter'],r.get('selection_maxiter','NA')])+' |')
 if stage=='development':
  lines += ['','## Bandwidth sensitivity (p=40; full MAVE)','','| Design | n | Bandwidth factor | Mean loss | Paired difference (MCSE) | Maxiter |','| --- | ---: | ---: | ---: | ---: | ---: |']
  for r in rows:
   if int(r['p'])==40 and r['method']=='full_mave':
    lines.append('| '+' | '.join([r['design'],r['n'],r['bandwidth_factor'],f"{float(r['mean_loss']):.5f}",f"{float(r['paired_delta']):+.5f} ({float(r['paired_mcse']):.5f})",r['maxiter']])+' |')
# Summarize actual iteration trajectories independently of final-fit summaries.
trajectory=defaultdict(list)
for path in (root/'confirmation').glob('trajectory-*.csv'):
 for r in csv.DictReader(path.open()):
  if int(r['p'])==40 and int(r['iteration']) in [1,5,10,20,40,60,100,150,200]:
   trajectory[r['design'],int(r['n']),int(r['iteration'])].append(float(r['loss']))
with (out/'trajectory.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['design','n','iteration','remaining_fits','mean_loss'])
 for (design,n,it),v in sorted(trajectory.items()):w.writerow([design,n,it,len(v),statistics.fmean(v)])
# Later-iteration means condition on not having stopped: do not treat as paired
# mean trajectories. Every final primary comparison includes all replications.
lines += ['','Trajectory tables condition on reaching the listed iteration. Compare them with final all-replicate summaries; they are not unconditional paired means.','']
(out/'comparisons.md').write_text('\n'.join(lines))
print('\n'.join(lines[:55]))
