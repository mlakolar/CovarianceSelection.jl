"""Validate support replay and combine its compact evaluation lookup."""
import csv,json,math,sys,tomllib
from pathlib import Path

planfile,out=sys.argv[1:];out=Path(out)
with open(planfile,'rb') as f:plan=tomllib.load(f)
count=4 if plan.get('candidate') else 22
expected=sum(s['replicates']*count for s in plan['scenario'])
seen=set();rows=[];maxerr=0.;maxleak=0.;checked=0
for path in sorted(out.glob('part-*.csv')):
    for r in csv.DictReader(path.open()):
        key=(r['scenario'],r['replicate'],r['id']);assert key not in seen;seen.add(key)
        assert int(r['extra'])==int(r['k'])-4+int(r['missing'])
        if r['qr_checked']=='true':
            checked+=1;maxerr=max(maxerr,float(r['qr_error']));maxleak=max(maxleak,float(r['direction_leak']))
        rows.append(r)
assert len(rows)==expected,(len(rows),expected)
with (out/'evaluation.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
receipt=dict(rows=len(rows),expected_rows=expected,duplicate_identities=0,qr_comparisons=checked,
             maximum_relative_ols_difference=maxerr,maximum_direction_projection_change=maxleak,
             support_cardinality_identity_verified=True)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
