"""Ensure the diagnostic correction leaves every primary reported quantity unchanged."""
import csv,json,math,sys
from pathlib import Path

root,parent=map(Path,sys.argv[1:]);checks=[]
for stage in ('development','confirmation_final'):
    def table(path):
        with path.open() as f:return {(r['scenario'],r['selector'],r['id']):r for r in csv.DictReader(f)}
    before=table(root/f'before-{stage}.csv');after=table(parent/stage/'summary.csv')
    assert before.keys()==after.keys()
    for key,old in before.items():
        new=after[key]
        for col in old:
            if col not in ('mean_missing','mean_extra'):assert old[col]==new[col],(stage,key,col)
        if new['selector']=='fixed' and '-screen-' in new['id']:
            k=int(new['id'].split('-')[-1])
            assert math.isclose(float(new['mean_extra']),k-4+float(new['mean_missing']),abs_tol=1e-12)
    checks.append(dict(stage=stage,groups=len(after),primary_quantities_unchanged=True,
                       selected_support_cardinality_identity=True))
(root/'correction_identity.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks))
