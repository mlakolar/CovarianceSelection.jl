"""Summarize fixed-penalty/bandwidth surfaces from saved pointwise results."""
import csv,json,math,statistics,sys
from collections import defaultdict
from pathlib import Path
root=Path(sys.argv[1]);rows=list(csv.DictReader((root/'summary.csv').open()))
groups=defaultdict(list)
for r in rows:
    if r['kind']=='pointwise':
        groups[tuple(r[k] for k in ['scenario','index_method','C','h_factor','pilot_factor','variant'])].append(r)
results=[]
for key,rr in sorted(groups.items()):
    out=dict(zip(['scenario','index_method','C','h_factor','pilot_factor','variant'],key))
    def avg(k):
        a=[float(r[k]) for r in rr if math.isfinite(float(r[k]))]
        return statistics.mean(a) if a else float('nan')
    def cover(a,b):return [float(r['coverage']) for r in rr if int(r['entry_a'])==a and int(r['entry_b'])==b]
    out.update(minimum_coverage=min(float(r['coverage']) for r in rr),mean_coverage=avg('coverage'),
        minimum_diagonal=min(cover(1,1)),maximum_diagonal=max(cover(1,1)),
        minimum_edge=min(cover(1,2)),maximum_edge=max(cover(1,2)),
        minimum_null=min(cover(1,4)),maximum_null=max(cover(1,4)),mean_width=avg('width'),
        failures=sum(int(r['failures']) for r in rr),index_exact=avg('index_exact'),
        nuisance_contains=avg('nuisance_contains'),nuisance_extra=avg('nuisance_extra'),
        mean_effective_n=avg('effective_n'),mean_lambda=avg('lambda_mean'))
    results.append(out)
with (root/'pointwise_surface.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(results[0]),lineterminator="\n");w.writeheader();w.writerows(results)
canonical=[r for r in results if float(r['C'])==1.1 and float(r['h_factor'])==1. and float(r['pilot_factor'])==2.]
(root/'canonical.json').write_text(json.dumps(canonical,indent=2)+'\n')
print(json.dumps(canonical,indent=2))
