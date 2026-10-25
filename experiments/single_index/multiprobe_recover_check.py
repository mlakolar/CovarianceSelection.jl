"""Test interrupted-record recovery, including retention of numerical failures."""
import csv,json,tempfile
from pathlib import Path
from multiprobe_recover import prepare,merge,write

root=Path(tempfile.mkdtemp(prefix='multiprobe-recovery-check-'))
plan=root/'plan.toml';plan.write_text('''seed = 7
shard_size = 2
validate = true
[[candidate]]
id = "test"
[[scenario]]
id = 1
n = 100
p = 40
design = "chain"
replicates = 2
''')
out=root/'out';out.mkdir(); recovery=root/'recovery'
fields=['scenario','replicate','n','p','design','selector','selected_id','status','loss']
def rows(rep):
    return [dict(scenario='1',replicate=str(rep),n='100',p='40',design='chain',selector=s,
                 selected_id='test',status='ok',loss='0.25') for s in ('fixed','validation','oracle')]
first=rows(1);first[0].update(status='failed',loss='NaN')
write(out/'part-0001.csv.tmp',first+rows(2)[:1],fields)
with (out/'part-0001.csv.tmp').open('a') as f:f.write('1,') # truncated final record
prepare(plan,out,recovery,8)
receipt=json.loads((recovery/'preparation.json').read_text())
assert receipt['retained_rows']==3 and receipt['missing_datasets']==1 and receipt['workers']==1
(recovery/'generated').mkdir();write(recovery/'generated/part-0002.csv',rows(2),fields)
merge(plan,out,recovery)
with (out/'part-0001.csv').open() as f: result=list(csv.DictReader(f))
assert len(result)==6 and sum(r['status']=='failed' for r in result)==1
before=(out/'part-0001.csv').read_bytes()
prepare(plan,out,root/'restart',8);merge(plan,out,root/'restart')
assert (out/'part-0001.csv').read_bytes()==before
assert json.loads((root/'restart/preparation.json').read_text())['missing_datasets']==0
print('Passed recovery checks: complete failed outcomes retained; incomplete datasets repeated; finalized shards unchanged.')
