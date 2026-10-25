"""Small numerical interpretation of validated tables, run under Slurm."""
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1]);sel=json.loads((root/'grid/selection.json').read_text())
rows=list(csv.DictReader((root/'reporting/summary.csv').open()));pairs=list(csv.DictReader((root/'reporting/paired.csv').open()))
roles={'opposed_scale':('diagonal_plus','edge'),'noncommuting':('diagonal','varying_edge'),'correlation':('diagonal_1','edge')}
fields=('coverage','availability','bias','empirical_sd','rms_se','length','index_loss','orientation_flip','missing','extra','effective_n','selected_mean')
result={'confirmation':{},'corrections':{},'sensitivities':{}}
for family,targets in roles.items():
 cid=sel['designs'][family]['development_best'];result['confirmation'][family]={}
 for n in (200,400,800,1600):
  cells=[r for r in rows if r['stage']=='confirmation' and r['design']==family and int(r['config'])==cid and int(r['n'])==n]
  result['confirmation'][family][n]=[{**{k:r[k] for k in ('point','target','index')},**{k:float(r[k]) for k in fields}} for r in cells]
  for target in targets:
   for index in ('known','estimated'):
    r=next(r for r in cells if float(r['point'])==.5 and r['target']==target and r['index']==index)
    print('MAIN',family,n,target,index,'cov',r['coverage'],'bias',r['bias'],'sd',r['empirical_sd'],'rmsse',r['rms_se'],'loss',r['index_loss'],'flip',r['orientation_flip'])
  pp=[r for r in pairs if r['stage']=='confirmation' and r['design']==family and int(r['n'])==n and int(r['config'])==cid and float(r['point'])==.5 and r['target']==targets[1]]
  print('PAIRED_EDGE',family,n,pp)
 corr=[r for r in rows if r['stage']=='corrections' and r['design']==family and float(r['point'])==0 and r['target']==targets[0]]
 result['corrections'][family]=[{**{k:r[k] for k in ('n','adjustment','target')},**{k:float(r[k]) for k in fields}} for r in corr]
 for n in (200,400,800,1600):print('CORRECTION_DIAGONAL_ZERO',family,n,[(r['adjustment'],r['coverage'],r['bias'],r['length']) for r in corr if int(r['n'])==n])
 sensitivity=[r for r in rows if r['stage'] in ('paired','double_training') and r['design']==family and int(r['config'])==cid and float(r['point'])==.5 and r['target']==targets[1]]
 result['sensitivities'][family]=[{**{k:r[k] for k in ('stage','n','ntrain','index')},**{k:float(r[k]) for k in fields}} for r in sensitivity]
 for r in sensitivity:print('SENSITIVITY_EDGE_HALF',family,r['stage'],r['n'],r['ntrain'],r['index'],r['coverage'],'loss',r['index_loss'],'flip',r['orientation_flip'])
(root/'reporting/interpretation.json').write_text(json.dumps(result,indent=2)+'\n')
