"""Pure plan generation. No local numerical simulations."""
from pathlib import Path
base=Path(__file__).parent
for stage,reps,ns,ps,hs,th,la,holdout in [
 ('smoke',2,[200,1600],[20,40],[1.], [1.],[.1],True),
 ('development',40,[200,400,800,1600],[20,40],[.7,1.,1.5],[.5,1.,1.5],[.03,.1,.25],True),
 ('confirmation',300,[200,400,800,1600],[20,40],[1.],[1.],[.1],True)]:
    rows=[f'seed = {826501 if stage=="development" else 2937511 if stage=="confirmation" else 713081}',
          f'shard_size = {2 if stage=="smoke" else 2}',
          'iterative_selection = true','max_updates = 60','dense_updates = 200','tolerance = 0.0001',f'bandwidth_factors = {hs}',
          f'threshold_levels = {th}',f'lasso_levels = {la}',
          f'holdout = {str(holdout).lower()}','holdout_bandwidth = 1.0']
    i=0
    for design in ['correlation','opposed_scale','fixed_diagonal','noncommuting']:
        for p in ps:
            for n in ns:
                i+=1
                bank='mixed' if design=='opposed_scale' else 'complete'
                rows += ['','[[scenario]]',f'id = {i}',f'design = "{design}"',
                         f'amplitude = {0.55 if design=="correlation" else 0.6}',
                         f'p = {p}',f'n = {n}',f'bank = "{bank}"',f'replicates = {reps}']
    (base/f'mave_{stage}.toml').write_text('\n'.join(rows)+'\n')
