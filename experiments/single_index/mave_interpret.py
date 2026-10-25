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
# Hybrid coefficients use the same feasible MAVE selection stage, not truth.
hybrid=root/'hybrid'/'summary.csv'
if hybrid.exists():
 hr=list(csv.DictReader(hybrid.open()));hg=defaultdict(dict)
 for r in hr:
  if r['method']=='stein_threshold_rrr_refit' and float(r['level'])==1.5:
   hg[r['design'],int(r['p']),int(r['n'])]['stein_threshold_rrr_refit_c15']=r
  if r['method']!='stein_threshold_rrr_refit' or float(r['level'])==1.:
   hg[r['design'],int(r['p']),int(r['n'])][r['method']]=r
 lines += ['','## Independent hybrid confirmation','', 'Both RRR and OLS refit probe slopes on the data-selected covariate support. The true-support RRR column is an infeasible diagnostic. These rows use the same 300 confirmation replications as the main comparison.','', '| Design | p | n | Stein10 | Threshold/RRR C=1 | Threshold/RRR C=1.5 | One sparse update/RRR | MAVE selection/RRR | MAVE selection/OLS | True-support RRR |','| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
 methods=['stein_screen10','stein_threshold_rrr_refit','stein_threshold_rrr_refit_c15','one_adaptive_select_rrr_refit','mave_select_rrr_refit','mave_select_ols_refit','oracle_support_rrr']
 for (design,p,n),g in sorted(hg.items()):
  lines.append('| '+' | '.join([design,str(p),str(n)]+[f"{float(g[m]['mean_loss']):.5f}" for m in methods])+' |')
 primary={};baseline={}
 for path in (root/'confirmation').glob('part-*.csv'):
  for r in csv.DictReader(path.open()):
   key=r['scenario'],r['replicate']
   if r['method']=='iterative_adaptive_lasso_refit':primary[key]=r
   if r['method']=='stein_screen10':baseline[key]=r
 hybrids={};checks=0
 for path in (root/'hybrid').glob('part-*.csv'):
  for r in csv.DictReader(path.open()):
   key=r['scenario'],r['replicate'];method=r['method']
   if method=='stein_screen10':
    assert abs(float(r['loss'])-float(baseline[key]['loss']))<1e-10;checks+=1
   if method=='mave_select_rrr_refit' and r['status']!='failed' and primary[key]['status']!='failed':
    q=primary[key]
    assert (r['missing'],r['extra'],r['selection_status'],r['selection_iterations'])==(q['missing'],q['extra'],q['selection_status'],q['selection_iterations']),(key,r,q)
    checks+=1
   if method=='stein_threshold_rrr_refit' and float(r['level'])==1.5:
    hybrids[key,'stein_threshold_rrr_refit_c15']=r
   if method!='stein_threshold_rrr_refit' or float(r['level'])==1.:
    hybrids[key,method]=r
 contrasts=defaultdict(list)
 for key,q in primary.items():
  a=hybrids[key,'mave_select_rrr_refit'];b=hybrids[key,'mave_select_ols_refit'];t=hybrids[key,'stein_threshold_rrr_refit']
  cell=a['design'],a['p'],a['n']
  contrasts[(*cell,'RRR hybrid minus MAVE coefficient refit')].append(float(a['loss'])-float(q['loss']))
  contrasts[(*cell,'RRR hybrid minus threshold/RRR')].append(float(a['loss'])-float(t['loss']))
  t15=hybrids[key,'stein_threshold_rrr_refit_c15']
  contrasts[(*cell,'RRR hybrid minus threshold/RRR C1.5')].append(float(a['loss'])-float(t15['loss']))
  contrasts[(*cell,'OLS hybrid minus RRR hybrid')].append(float(b['loss'])-float(a['loss']))
 with (out/'hybrid-paired.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['design','p','n','contrast','replicates','paired_delta','paired_mcse'])
  for k,v in sorted(contrasts.items()):w.writerow([*k,len(v),statistics.fmean(v),statistics.stdev(v)/math.sqrt(len(v))])
 (out/'hybrid-validation.json').write_text(json.dumps(dict(paired_checks=checks,primary_cases=len(primary),hybrid_summary_rows=len(hr),checks='Baseline losses agree; successful MAVE-selected support counts and selection-stage status/iteration counts agree across coefficient-refit arms; all planned pairs retained.'),indent=2)+'\n')

# Validate stored trajectories against every successful final dense fit.
trajectory_checks={}
for stage in ['development','confirmation']:
 expected={}
 for path in (root/stage).glob('part-*.csv'):
  for r in csv.DictReader(path.open()):
   if r['method']=='full_mave' and r['status']!='failed':
    expected[r['scenario'],r['replicate'],r['bandwidth_factor']]=(int(r['iterations']),float(r['loss']))
 observed={};records=0
 for path in (root/stage).glob('trajectory-*.csv'):
  for r in csv.DictReader(path.open()):
   key=r['scenario'],r['replicate'],r['bandwidth_factor'];it=int(r['iteration']);loss=float(r['loss'])
   assert key in expected and it==observed.get(key,(0,0.))[0]+1,(stage,key,it)
   assert math.isfinite(loss) and -1e-9<=loss<=1+1e-9
   observed[key]=(it,loss);records+=1
 assert set(expected)==set(observed),(stage,len(expected),len(observed))
 for key,(it,loss) in expected.items():
  assert observed[key][0]==it and abs(observed[key][1]-loss)<1e-10,(stage,key)
 trajectory_checks[stage]=dict(successful_dense_fits=len(expected),trajectory_records=records,checks='Contiguous iteration numbering, unique fit/iteration keys, finite bounded losses, exact final-iteration count and agreement with primary final direction loss.')
(out/'trajectory-validation.json').write_text(json.dumps(trajectory_checks,indent=2)+'\n')

diagpath=root/'convergence'/'summary.csv'
if diagpath.exists():
 diag=list(csv.DictReader(diagpath.open()))
 lines += ['','## Bounded longer-convergence diagnostic','',f"All {len(diag)} cases used the first two development replications for each design at p=40,n=400/800, fixed before inspecting their errors. {sum(r['status']=='converged' for r in diag)} converged within 2000 updates. {sum(float(r['loss2000'])>float(r['stein_loss']) for r in diag)} ended with larger direction loss than their Stein initializers. This small diagnostic does not replace the all-replicate primary comparison.",'']
 (out/'convergence-diagnostic.json').write_text(json.dumps(dict(cases=len(diag),converged=sum(r['status']=='converged' for r in diag),worse_than_stein=sum(float(r['loss2000'])>float(r['stein_loss']) for r in diag),maximum_iterations=max(int(r['iterations']) for r in diag)),indent=2)+'\n')

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
import subprocess,os
plot_python=os.environ.get("SIG_PLOT_PYTHON","/project2/mkolar_1314/SingleIndexGraph/20261024-sig-finalvalidation/.report-env/bin/python")
subprocess.run([plot_python,str(root/'experiments/single_index/mave_figures.py'),str(root)],check=True)
