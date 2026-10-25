"""Freeze specifications using known-beta development only, before confirmation."""
import csv,json,math,sys,tomllib
from pathlib import Path
from collections import defaultdict
from promising_coverage_plans import write,scenarios

out=Path(sys.argv[1]);plan=tomllib.loads(Path('experiments/single_index/promising_coverage_grid.toml').read_text())
rows=list(csv.DictReader((out/'summary.csv').open()))
if any(r['index']!='known' or r['adjustment']!='basic' for r in rows):raise ValueError('selection requires basic known-beta results only')
configs={c['id']:c for c in plan['candidate']};groups=defaultdict(list)
for r in rows:groups[r['design'],int(r['config'])].append(r)
ranked=defaultdict(list)
for (design,cid),rr in groups.items():
    coverage=[float(r['coverage']) for r in rr]
    failures=sum(1-float(r['availability']) for r in rr)/len(rr)
    # All entry/point/sample-size cells count equally. Penalize failures explicitly.
    score=sum(abs(c-.95) for c in coverage)/len(rr)+2*failures
    ranked[design].append((score,cid))
selected=[];receipt={};poor=[]
for design,scores in sorted(ranked.items()):
    scores.sort();best=scores[0][1]
    baseline=next(c['id'] for c in plan['candidate'] if c['C']==1.1 and c['hfactor']==1.5 and c['pilot_ratio']==2.)
    ids=list(dict.fromkeys([baseline,best]))
    if len(ids)==1:ids.append(next(cid for score,cid in scores if cid!=best))
    for cid in ids:selected.append(dict(configs[cid],design=design))
    chosen=groups[design,best]
    bad=[r for r in chosen if int(r['n'])>=400 and (float(r['coverage'])<.925 or float(r['availability'])<.99)]
    if bad:poor.append(design)
    receipt[design]={'baseline':baseline,'development_best':best,'selected':ids,'ranking':scores,'poor_known_beta_cells':[{k:r[k] for k in ('n','point','target','coverage','availability')} for r in bad]}
# Seed namespaces deliberately separated from all development and canary seeds.
write('paired',scenarios(300),selected,['known','estimated','estimated_full'],121736511,10)
write('confirmation',scenarios(2000),selected,['known','estimated'],221736511,20)
# Additional index training observations, without changing precision n.
write('double_training',scenarios(300,(400,800),2),selected,['known','estimated'],321736511,10)
corrections=[]
for design in poor:
    cid=receipt[design]['development_best']
    for adjustment in ('basic','df','hc3','inverse'):
        corrections.append(dict(configs[cid],id=cid+{'basic':0,'df':100,'hc3':200,'inverse':300}[adjustment],design=design,adjustment=adjustment))
if corrections:
    write('corrections',[s for s in scenarios(300) if s['design'] in poor],corrections,['known'],421736511,10)
selection={'criterion':'Mean absolute known-beta coverage error across all entry/point/n cells plus twice failure rate; fixed baseline retained. No estimated-beta results used.','poor_threshold':'On development best: coverage <0.925 or availability <0.99 at n>=400, any requested point/entry. Exploratory trigger, not a significance test.','designs':receipt,'corrections_designs':poor}
(out/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
print(json.dumps(selection,indent=2))
