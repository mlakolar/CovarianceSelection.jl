"""Paired first-stage diagnostics, retaining empty-support and anchor failures."""
import csv,json,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path

root=Path(sys.argv[1]); groups=defaultdict(list);seen=set()
for r in csv.DictReader((root/'comparison.csv').open()):
    key=(r['n'],r['zero_first'],r['policy']);identity=key+(r['rep'],)
    assert identity not in seen
    seen.add(identity);groups[key].append(r)
assert len(seen)==4000
rows=[]
for key,rr in sorted(groups.items()):
    assert len(rr)==500
    failures=Counter(r['failure_kind'] for r in rr)
    ok=[r for r in rr if not r['failure_kind']]
    rows.append(dict(zip(('n','zero_first','policy'),key))|dict(replicates=500,
        available=len(ok),empty_support=failures['empty_support'],missing_anchor=failures['missing_anchor'],
        other_failures=sum(bool(r['failure_kind']) and r['failure_kind'] not in ('empty_support','missing_anchor') for r in rr),
        error2=statistics.mean(float(r['error2']) for r in ok) if ok else float('nan'),
        projective_error2=statistics.mean(float(r['projective_error2']) for r in ok) if ok else float('nan')))
with (root/'summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
receipt=dict(rows=len(seen),expected_rows=4000,groups=8,replicates_per_group=500,duplicate_identities=0,failures_retained=True)
(root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(receipt=receipt,summary=rows),indent=2))
