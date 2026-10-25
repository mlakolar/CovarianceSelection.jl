"""Coverage ranges across points, with dimensions and failures kept explicit."""
import csv,json,sys
from collections import defaultdict
from pathlib import Path

def report(root):
    root=Path(root); rr=list(csv.DictReader((root/'summary.csv').open()))
    groups=defaultdict(list)
    for r in rr:
        groups[tuple(r[k] for k in ('scenario','index_method','variant','entry_a','entry_b'))].append(r)
    ranges=[]
    for key,rows in sorted(groups.items()):
        out=dict(zip(('scenario','index_method','variant','entry_a','entry_b'),key))
        out.update(coverage_min=min(float(r['coverage']) for r in rows),
            coverage_max=max(float(r['coverage']) for r in rows),
            maximum_failures=max(int(r['failures']) for r in rows),
            final_effective_n_min=min(float(r['final_effective_n']) for r in rows),
            final_effective_n_max=max(float(r['final_effective_n']) for r in rows),
            nuisance_contains_min=min(float(r['nuisance_contains']) for r in rows))
        ranges.append(out)
    with (root/'coverage_ranges.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(ranges[0]),lineterminator='\n');w.writeheader();w.writerows(ranges)
    primary=[r for r in ranges if r['variant']=='union_refit']
    for r in primary:
        print(r['scenario'],r['index_method'],(r['entry_a'],r['entry_b']),
            f"coverage={r['coverage_min']:.3f}--{r['coverage_max']:.3f}",
            f"failures<={r['maximum_failures']}",f"n_eff={r['final_effective_n_min']:.1f}--{r['final_effective_n_max']:.1f}")

if __name__=='__main__': report(sys.argv[1])
