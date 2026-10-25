"""Known directions must give identical fixed/full bootstrap endpoints."""
import csv,math,sys
from pathlib import Path

rows=[]
for p in Path(sys.argv[1]).glob('shard-*.csv'):
    rows.extend(csv.DictReader(p.open()))
fixed={(r['scenario'],r['rep'],r['entry_a'],r['entry_b'],r['point']):r
    for r in rows if r['index_method']=='known' and r['interval_method']=='fixed_index'}
full={(r['scenario'],r['rep'],r['entry_a'],r['entry_b'],r['point']):r
    for r in rows if r['index_method']=='known' and r['interval_method']=='full_pipeline'}
assert fixed.keys()==full.keys() and len(fixed)==18
for key,r in fixed.items():
    s=full[key]
    for field in ('lower','upper','covered','failure','bootstrap_failures'):
        assert r[field]==s[field],(key,field,r[field],s[field])
assert any(not r['failure'] for r in fixed.values())
print('Passed paired-bootstrap identity checks at 18 entry/location/dataset combinations.')
