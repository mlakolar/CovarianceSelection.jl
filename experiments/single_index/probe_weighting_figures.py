"""Plot fixed, prespecified comparisons. Run only under Discovery Slurm."""
import csv,json,sys,math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);out=root/'reporting';out.mkdir(exist_ok=True)
rows=list(csv.DictReader((root/'confirmation/summary.csv').open()))
main=[r for r in rows if float(r['signal'])==1 and r['beta_kind']=='equal']
styles=[('ols',20,1.,'Current OLS, all 20'),('rrr',20,1.,'RRR, all 20'),('weighted_rrr',20,.2,'Weighted RRR, all 20'),('weighted_rrr',10,.2,'Weighted RRR, select 10'),('learned_crossfit',20,.2,'Learned probe, cross-fit')]
for bank in ('projections','mixed','complete'):
    fig,axes=plt.subplots(2,3,figsize=(13,7),sharex=True,sharey=True)
    for i,design in enumerate(('chain','pairs')):
        for j,p in enumerate((20,40,160)):
            ax=axes[i,j]
            if bank=='complete' and p==160:ax.axis('off');continue
            for method,k,rho,label in styles:
                a=sorted([r for r in main if r['design']==design and int(r['p'])==p and r['bank']==bank and r['method']==method and int(r['k'])==k and float(r['rho'])==rho],key=lambda r:int(r['n']))
                ax.errorbar([int(r['n']) for r in a],[float(r['mean_loss']) for r in a],yerr=[1.96*float(r['mcse_loss']) for r in a],label=label,marker='o',ms=3,capsize=2)
            for method,label,color in [('ols','Known informative probe','gray'),('true_support_ols','Known probe + true support','black')]:
                a=sorted([r for r in main if r['design']==design and int(r['p'])==p and r['bank']=='oracle_probe' and r['method']==method],key=lambda r:int(r['n']))
                ax.plot([int(r['n']) for r in a],[float(r['mean_loss']) for r in a],ls='--',color=color,label=label)
            ax.set_title(f'{design}, p={p}');ax.set_xscale('log',base=2);ax.set_ylim(0,1);ax.grid(alpha=.2)
            if j==0:ax.set_ylabel('Mean squared-sine error')
            if i==1:ax.set_xlabel('Sample size')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
    fig.suptitle(f'{bank}: fixed aggregation rules; fresh confirmation')
    fig.tight_layout(rect=(0,.12,1,.96))
    for ext in ('png','pdf'):fig.savefig(out/f'{bank}_direction.{ext}',dpi=160)
    plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True)
for ax,design in zip(axes,('chain','pairs')):
    for bank in ('projections','mixed','complete'):
        for method,k,rho,label in (styles[0],styles[2],styles[4]):
            a=sorted([r for r in main if r['design']==design and int(r['p'])==40 and r['bank']==bank and r['method']==method and int(r['k'])==k and float(r['rho'])==rho],key=lambda r:int(r['n']))
            ax.plot([int(r['n']) for r in a],[float(r['mean_weight_l1']) for r in a],label=f'{bank}: {label}',marker='o',ms=3)
    ax.set_title(f'{design}, p=40');ax.set_xscale('log',base=2);ax.set_xlabel('Sample size');ax.grid(alpha=.2)
axes[0].set_ylabel('Mean L1 distance between kernel weights')
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=7)
fig.suptitle('Neighborhood distortion at h = n^(-1/5); no precision coverage claim')
fig.tight_layout(rect=(0,.2,1,.92))
for ext in ('png','pdf'):fig.savefig(out/f'neighborhood_distortion.{ext}',dpi=160)
plt.close(fig)
import shutil
for stage in ('canary','development','confirmation'):
    dest=out/stage;dest.mkdir(exist_ok=True)
    for name in ('summary.csv','paired.csv','receipt.json','failures.csv'):
        path=root/stage/name
        if path.exists():shutil.copy2(path,dest/name)
    for path in (root/stage).glob('inputs-*.sha256'):shutil.copy2(path,dest/path.name)
# Outcome report separates controls and feasible methods; no selection rule.
text=['# Fixed aggregation comparison','', 'All fixed rules were evaluated on fresh confirmation datasets.','No precision intervals were constructed. Full results and paired differences are in the CSV files.','', '| Design | p | n | Bank | Method | k | rho | Mean loss | MCSE | Failures |','|---|---:|---:|---|---|---:|---:|---:|---:|---:|']
for r in sorted(main,key=lambda r:(r['design'],int(r['p']),r['bank'],r['method'],int(r['k']),float(r['rho']))):
    if int(r['n'])!=3200:continue
    if r['bank']!='oracle_probe' and (r['method'],int(r['k']),float(r['rho'])) not in [x[:3] for x in styles]:continue
    text.append('| '+' | '.join([r['design'],r['p'],r['n'],r['bank'],r['method'],r['k'],r['rho'],f"{float(r['mean_loss']):.4f}",f"{float(r['mcse_loss']):.4f}",r['failures']])+' |')
(out/'RESULTS.md').write_text('\n'.join(text)+'\n')
print(json.dumps({'figures':8,'reporting':str(out)}))
