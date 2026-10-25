"""Recover complete replicate records from interrupted atomic shards.

Only I/O batching changes. The original fitter, scenario ids, seeds, and probe
seeds are reused. A numerical failure is a completed outcome, never a retry.
"""
import csv,json,math,sys,tomllib
from collections import defaultdict
from pathlib import Path

def load(planfile):
    with open(planfile,'rb') as f:return tomllib.load(f)

def tasks(plan):
    return [(s,range(first,min(first+plan['shard_size'],s['replicates']+1)))
            for s in plan['scenario'] for first in range(1,s['replicates']+1,plan['shard_size'])]

def labels(plan):
    return {('fixed',c['id']) for c in plan['candidate']}|{('validation','validation'),('oracle','oracle')}

def key(row):
    return row['selector'],row['selected_id'] if row['selector']=='fixed' else row['selector']

def write(path,rows,fields):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)

def valid_group(rows,expected,s,rep,fields):
    if len(rows)!=len(expected) or {key(r) for r in rows}!=expected:return False
    for r in rows:
        if set(r)!=set(fields) or any(v is None for v in r.values()):return False
        if r['scenario']!=str(s['id']) or r['replicate']!=str(rep):return False
        if r['n']!=str(s['n']) or r['p']!=str(s['p']) or r['design']!=s['design']:return False
        if r['status'] not in ('ok','failed'):return False
        if r['status']=='ok':
            try:x=float(r['loss'])
            except ValueError:return False
            if not math.isfinite(x) or not 0<=x<=1:return False
    return True

def prepare(planfile,out,root,workers):
    out=Path(out);root=Path(root);root.mkdir(exist_ok=True)
    plan=load(planfile); expected=labels(plan);work=[]; retained=[]; metadata=[]; fields=None
    # Each original scenario occurs once and uses its original id. One-rep
    # task indices are offsets in the UNCHANGED full scenario sequence.
    offsets={};offset=0
    for s in plan['scenario']:offsets[s['id']]=offset;offset+=s['replicates']
    for task,(s,reps) in enumerate(tasks(plan),1):
        final=out/f'part-{task:04}.csv'
        if final.exists():continue
        tmp=Path(str(final)+'.tmp'); groups=defaultdict(list)
        if tmp.exists():
            with tmp.open() as f:
                reader=csv.DictReader(f);fields=reader.fieldnames
                for row in reader:
                    if None not in row and (row.get('replicate') or '').isdigit():
                        groups[int(row['replicate'])].append(row)
        keep=[];missing=[]
        for rep in reps:
            rr=groups.get(rep,[])
            if fields and valid_group(rr,expected,s,rep,fields):keep.extend(rr)
            else:
                one=offsets[s['id']]+rep;work.append(one);missing.append((rep,one))
        if fields:write(root/f'retained-{task:04}.csv',keep,fields)
        retained.extend(keep)
        metadata.append(dict(task=task,scenario=s['id'],reps=list(reps),missing=missing))
    original=Path(planfile).read_text()
    # The generated plans place shard_size on its own line exactly once.
    lines=original.splitlines();count=sum(l.startswith('shard_size =') for l in lines)
    assert count==1
    (root/'one_rep.toml').write_text('\n'.join('shard_size = 1' if l.startswith('shard_size =') else l for l in lines)+'\n')
    one=load(root/'one_rep.toml');assert one.pop('shard_size')==1
    base=plan.copy();base.pop('shard_size');assert one==base
    workers=min(workers,len(work))
    for k in range(workers):
        (root/f'worker-{k+1}.txt').write_text(''.join(f'{i}\n' for i in work[k::workers]))
    (root/'workers.txt').write_text(''.join(f'{k}\n' for k in range(1,workers+1)))
    receipt=dict(original_plan=str(planfile),scientific_plan_identical=True,
                 batching_only_change=f"shard_size {plan['shard_size']} -> 1",unfinished_shards=metadata,
                 retained_rows=len(retained),missing_datasets=len(work),workers=workers,fields=fields)
    (root/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('retained_rows','missing_datasets','workers')}))

def merge(planfile,out,root):
    out=Path(out);root=Path(root);plan=load(planfile);expected=labels(plan)
    state=json.loads((root/'preparation.json').read_text()); fields=state['fields']; merged=0
    for item in state['unfinished_shards']:
        task=item['task'];s=next(s for s in plan['scenario'] if s['id']==item['scenario'])
        keep=root/f'retained-{task:04}.csv'
        rows=list(csv.DictReader(keep.open())) if keep.exists() else []
        for rep,one in item['missing']:
            path=root/'generated'/f'part-{one:04}.csv'
            with path.open() as f:
                reader=csv.DictReader(f)
                if fields is None:fields=reader.fieldnames
                added=list(reader)
            assert valid_group(added,expected,s,rep,fields),(task,rep)
            rows.extend(added)
        groups=defaultdict(list)
        for r in rows:groups[int(r['replicate'])].append(r)
        assert set(groups)==set(item['reps'])
        for rep,rr in groups.items():assert valid_group(rr,expected,s,rep,fields)
        final=out/f'part-{task:04}.csv';assert not final.exists()
        temporary=Path(str(final)+'.recovered')
        write(temporary,sorted(rows,key=lambda r:(int(r['replicate']),key(r))),fields)
        temporary.rename(final);merged+=1
    (root/'receipt.json').write_text(json.dumps(dict(merged_shards=merged,
        retained_rows=state['retained_rows'],recomputed_datasets=state['missing_datasets'],
        scientific_plan_identical=True,failed_outcomes_retained=True),indent=2)+'\n')
    print('Recovered shards:',merged)

if __name__=='__main__':
    action,plan,out,root,*rest=sys.argv[1:]
    if action=='prepare':prepare(plan,out,root,int(rest[0]))
    elif action=='merge':merge(plan,out,root)
    else:raise ValueError(action)
