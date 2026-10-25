"""Freeze independent sensitivity plans after the initial coverage ladder."""
import csv,json,sys,tomllib
from pathlib import Path
from gaussian_fixed_plans import save

def followup(planpath,summarypath,destination):
    destination=Path(destination)
    plan=tomllib.loads(Path(planpath).read_text())
    rows=list(csv.DictReader(Path(summarypath).open()))
    byname={c['name']:c for c in plan['scenario']}
    families={}
    for c in plan['scenario']:
        families.setdefault((c['graph'],c['p']),[]).append(c)
    selected=[]; sensitivity=[]
    for key,cases in sorted(families.items()):
        cases.sort(key=lambda c:c['n'])
        minima=[]
        for c in cases:
            rr=[r for r in rows if r['scenario']==c['name'] and r['variant']=='union_refit']
            assert len(rr)==9
            minima.append(min(float(r['coverage']) for r in rr))
        stable=[i for i in range(len(cases)) if all(x>=.93 for x in minima[i:])]
        i=stable[0] if stable else len(cases)-1
        indices=sorted(set((max(0,i-1),i)))
        if len(indices)==1:indices=[0,1]
        chosen=[cases[j] for j in indices]
        selected.append(dict(graph=key[0],p=key[1],coverage_minima=minima,
            selected_n=[c['n'] for c in chosen],stable_93_percent_observed=bool(stable)))
        for c in chosen:
            # Change only one bandwidth at a time, including a fresh baseline.
            for hf,pf in ((.7,2.),(1.,2.),(1.3,2.),(1.,1.),(1.,3.)):
                sensitivity.append(dict(c,id=len(sensitivity)+1,
                    name=f"{c['name']}_h{hf}_b{pf}",replicates=500,
                    h_factors=[hf],pilot_factors=[pf]))
    save(destination/'sample_size_bandwidth.toml',
        dict(stage='pointwise',seed_base=947000000,shard_size=50,penalty_probability=.1),sensitivity)
    (destination/'sample_size_followup_selection.json').write_text(json.dumps(selected,indent=2)+'\n')
    print(json.dumps(selected,indent=2))

if __name__=='__main__':followup(*sys.argv[1:])
