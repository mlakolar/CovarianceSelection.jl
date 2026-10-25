"""Availability, coverage ranges, and paired changes for bootstrap intervals."""
import csv,math,statistics,sys
from collections import defaultdict
from pathlib import Path

def report(root):
    root=Path(root);summary=list(csv.DictReader((root/'summary.csv').open()))
    groups=defaultdict(list)
    for r in summary:
        groups[(r['scenario'],r['index_method'],r['interval_method'])].append(r)
    ranges=[]
    for key,rr in sorted(groups.items()):
        out=dict(zip(('scenario','index_method','interval_method'),key))
        out.update(coverage_min=min(float(r['coverage']) for r in rr),
            coverage_max=max(float(r['coverage']) for r in rr),
            available_min=min((int(r['n'])-int(r['failures']))/int(r['n']) for r in rr),
            available_max=max((int(r['n'])-int(r['failures']))/int(r['n']) for r in rr))
        ranges.append(out);print(out)
    with (root/'coverage_ranges.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(ranges[0]),lineterminator='\n');w.writeheader();w.writerows(ranges)
    cells=defaultdict(dict)
    for path in root.glob('shard-*.csv'):
        for r in csv.DictReader(path.open()):
            key=tuple(r[k] for k in ('scenario','index_method','entry_a','entry_b','point','rep'))
            cells[key][r['interval_method']]=r
    paired=defaultdict(list)
    for key,rr in cells.items():
        assert set(rr)=={'wald','fixed_index','full_pipeline'}
        for mode in ('fixed_index','full_pipeline'):
            paired[key[:-1]+(mode,)].append(int(rr[mode]['covered'])-int(rr['wald']['covered']))
    with (root/'paired_coverage.csv').open('w') as f:
        w=csv.writer(f,lineterminator='\n')
        w.writerow(['scenario','index_method','entry_a','entry_b','point','interval_method','n','coverage_change','paired_mcse'])
        for key,dd in sorted(paired.items()):
            w.writerow([*key,len(dd),statistics.mean(dd),statistics.stdev(dd)/math.sqrt(len(dd))])

if __name__=='__main__':report(sys.argv[1])
