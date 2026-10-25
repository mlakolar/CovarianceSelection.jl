"""Truth-aligned sign diagnostic ONLY; never a feasible CI or primary result."""
import csv,json,math,sys
from collections import defaultdict
from pathlib import Path
root=Path(sys.argv[1]);sel=json.loads((root/'grid/selection.json').read_text());sums=defaultdict(lambda:defaultdict(float))
for path in sorted((root/'confirmation').glob('part-*.csv')):
 groups=defaultdict(dict)
 for r in csv.DictReader(path.open()):
  if r['index']!='estimated' or int(r['config'])!=sel['designs'][r['design']]['development_best']:continue
  groups[(r['scenario'],r['replicate'],r['config'],r['target'])][float(r['point'])]=r
 for ident,points in groups.items():
  if set(points)!={-.5,0.,.5}:raise ValueError(f'incomplete point grid {ident}')
  for point,r in points.items():
   flip=int(r['orientation_flip']);aligned=points[-point] if flip else r
   good=aligned['status']=='ok';truth=float(r['truth'])
   covered=good and abs(float(aligned['estimate'])-truth)<=1.959963984540054*float(aligned['se'])
   original=int(r['covered']);diff=int(covered)-original
   s=sums[r['design'],r['n'],r['target'],point];s['replicates']+=1;s['original']+=original;s['aligned']+=covered;s['flips']+=flip;s['flipped_covered']+=flip*original;s['unflipped_covered']+=(1-flip)*original;s['diff2']+=diff*diff
rows=[]
for (design,n,target,point),s in sorted(sums.items()):
 R=s['replicates'];F=s['flips'];diff=s['aligned']-s['original']
 if R!=2000:raise ValueError('wrong replication count')
 rows.append(dict(design=design,n=int(n),target=target,point=point,replicates=int(R),original_coverage=s['original']/R,truth_aligned_coverage=s['aligned']/R,orientation_flip_rate=F/R,coverage_given_flip=s['flipped_covered']/F if F else None,coverage_given_no_flip=s['unflipped_covered']/(R-F) if R>F else None,alignment_difference_mcse=math.sqrt(max(0,(s['diff2']-diff*diff/R)/(R-1)/R))))
receipt=dict(scope='Infeasible sign diagnostic only. Flip means estimated direction has negative inner product with true beta. Swap the existing +/- point fits; no refitting and no change to primary results.',rows=rows)
(root/'reporting/orientation_diagnostic.json').write_text(json.dumps(receipt,indent=2)+'\n')
for r in rows:
 if r['point']==.5 and r['target'] in ('edge','varying_edge'):print(json.dumps(r))
