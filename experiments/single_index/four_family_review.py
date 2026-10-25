"""Targeted checks of observed fit failures and reported support diagnostics."""
import csv,json,sys,math
from pathlib import Path
from collections import Counter
root=Path(sys.argv[1]);out=root/'reporting'
rows=list(csv.DictReader((root/'comparison/summary.csv').open()))
failed=list(csv.DictReader((root/'comparison/failures.csv').open())) if (root/'comparison/failures.csv').exists() else []
for r in failed:
    assert r['method']=='split',r
    assert 'no positive empirical signal eigenvalue' in r['reason'],r
    assert float(r['loss'])==1.,r
assert sum(int(r['failures']) for r in rows)==len(failed)
for r in rows:
    if math.isfinite(float(r['mean_extra'])):
        assert abs(float(r['mean_extra'])-(int(r['k'])-4+float(r['mean_missing'])))<1e-8,r
lookup={(r['scenario'],r['bank'],r['method'],r['k'],r['rho']):r for r in rows}
matched=0
for r in rows:
    if r['d']=='4' and r['bank']=='oracle_probe' and r['method']=='ols':
        b=lookup[(r['scenario'],'oracle_probe','true_support_ols','4','0.0')]
        assert abs(float(r['mean_loss'])-float(b['mean_loss']))<1e-10
        matched+=1
counts=Counter((r['design'],r['amplitude'],r['p'],r['d'],r['n'],r['bank'],r['method']) for r in failed)
with (out/'failure_summary.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['design','amplitude','p','d','n','bank','method','failed_fits'])
    for key,count in sorted(counts.items()):w.writerow([*key,count])
receipt=dict(failed_fits=len(failed),all_failures_split_moment_nonpositive_signal=True,all_failures_retained_with_maximum_loss=True,support_cardinality_identities=True,identical_d4_oracle_controls=matched,canary_comparison_input_hashes_identical=(root/'canary/inputs-13445394.sha256').read_text()==(root/'comparison/inputs-13445513.sha256').read_text())
assert receipt['canary_comparison_input_hashes_identical']
(out/'review_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
for r in rows:
    if r['p']=='40' and r['n'] in ('400','800') and r['d']=='20' and r['amplitude'] in ('0.55','0.6') and r['bank']=='complete' and r['method'] in ('rrr','weighted_rrr','learned_crossfit') and (r['method']!='weighted_rrr' or r['rho']=='0.2'):
        print(*(r[k] for k in ('design','n','method','k','mean_loss','mcse_loss','mean_bandwidth_ratio','mean_weight_l1')))
