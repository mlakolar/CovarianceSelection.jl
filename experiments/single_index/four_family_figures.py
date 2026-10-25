"""Final four-design reporting; runs under Discovery Slurm."""
import csv,json,sys,shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);dest=root/'reporting';dest.mkdir(exist_ok=True)
rows=list(csv.DictReader((root/'comparison/summary.csv').open()))
styles=[('complete','ols',1.,'full','Complete OLS'),('complete','rrr',1.,'full','Complete RRR'),('complete','weighted_rrr',.2,'full','Complete weighted RRR (.2)'),('complete','weighted_rrr',.2,'screen','Complete weighted RRR (.2), select 10'),('mixed','rrr',1.,'screen','Mixed RRR, at most 10'),('complete','learned_crossfit',.2,'full','Complete learned probe'),('oracle_probe','ols',0.,'full','Known informative probe'),('oracle_probe','true_support_ols',0.,'true','Known probe + true support')]
styles.append(('complete','rrr',1.,'screen','Complete RRR, select 10'))
for family in ('correlation','opposed_scale','fixed_diagonal','noncommuting'):
    amps=sorted({float(r['amplitude']) for r in rows if r['design']==family})
    fig,axes=plt.subplots(2,3,figsize=(13,7),sharex=True,sharey=True)
    for i,d in enumerate((4,20)):
        for j,a in enumerate(amps):
            ax=axes[i,j]
            for style_index,(bank,method,rho,kind,label) in enumerate(styles):
                if d==4 and kind=='screen' and bank=='complete':continue
                k=4 if kind=='true' else (10 if kind=='screen' and d==20 else d)
                subset=sorted([r for r in rows if r['design']==family and int(r['p'])==40 and int(r['d'])==d and float(r['amplitude'])==a and r['bank']==bank and r['method']==method and int(r['k'])==k and float(r['rho'])==rho],key=lambda r:int(r['n']))
                color='black' if kind=='true' else 'gray' if bank=='oracle_probe' else f'C{style_index}'
                ax.errorbar([int(r['n']) for r in subset],[float(r['mean_loss']) for r in subset],yerr=[1.96*float(r['mcse_loss']) for r in subset],label=label,marker='o',ms=3,capsize=2,ls='--' if bank=='oracle_probe' else '-',color=color)
            ax.set_title(f'd={d}, amplitude={a:g}');ax.set_xscale('log',base=2);ax.set_ylim(0,1);ax.grid(alpha=.2)
            if j==0:ax.set_ylabel('Mean squared-sine direction error')
            if i==1:ax.set_xlabel('Sample size')
    handles,labels=axes[1,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
    fig.suptitle(f'{family}: response dimension 40; 300 paired repetitions')
    fig.tight_layout(rect=(0,.15,1,.96))
    for ext in ('png','pdf'):fig.savefig(dest/f'{family}_direction.{ext}',dpi=160)
    plt.close(fig)
# The scalar baselines are feasible and use no graph information.
fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
for ax,family in zip(axes.ravel(),('correlation','opposed_scale','fixed_diagonal','noncommuting')):
    a=max(float(r['amplitude']) for r in rows if r['design']==family)
    for bank,method,rho,label in [('global_trace','ols',1.,'Trace probe'),('global_products','ols',1.,'All-products probe'),('complete','rrr',1.,'Complete RRR'),('complete','weighted_rrr',.2,'Complete weighted RRR (.2)'),('mixed','second',1.,'Mixed second-order Stein')]:
        subset=sorted([r for r in rows if r['design']==family and int(r['p'])==40 and int(r['d'])==20 and float(r['amplitude'])==a and r['bank']==bank and r['method']==method and int(r['k'])==20 and float(r['rho'])==rho],key=lambda r:int(r['n']))
        ax.plot([int(r['n']) for r in subset],[float(r['mean_loss']) for r in subset],marker='o',ms=3,label=label)
    ax.set_title(f'{family}, a={a:g}');ax.set_xscale('log',base=2);ax.set_ylim(0,1);ax.grid(alpha=.2)
    ax.set_xlabel('Sample size');ax.set_ylabel('Mean squared-sine error')
handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
fig.suptitle('Generic scalar probes and second-order comparison; d=20, p=40')
fig.tight_layout(rect=(0,.12,1,.94))
for ext in ('png','pdf'):fig.savefig(dest/f'scalar_and_second_order.{ext}',dpi=160)
plt.close(fig)
for stage in ('canary','comparison'):
    target=dest/stage;target.mkdir(exist_ok=True)
    for name in ('summary.csv','paired.csv','receipt.json','failures.csv'):
        p=root/stage/name
        if p.exists():shutil.copy2(p,target/name)
    for p in (root/stage).glob('inputs-*.sha256'):shutil.copy2(p,target/p.name)
shutil.copytree(root/'population',dest/'population',dirs_exist_ok=True)
shutil.copy2(root/'population-inputs.sha256',dest/'population-inputs.sha256')
text=['# Four-design index comparison','', 'All methods and amplitudes were prespecified. No precision estimation or inference ran.','The table shows strongest amplitude, response dimension 40, n=400 and n=800.','All settings, fixed shrinkage levels, tail errors, failures and paired controls are in CSV files.','', '| Family | a | d | n | Probe bank | Method | k | rho | Mean loss | MCSE |','|---|---:|---:|---:|---|---|---:|---:|---:|---:|']
for r in rows:
    if int(r['p'])!=40 or int(r['n']) not in (400,800):continue
    if float(r['amplitude'])!=max(float(x['amplitude']) for x in rows if x['design']==r['design']):continue
    if not (r['bank']=='oracle_probe' or (r['bank']=='complete' and r['method']=='rrr') or (r['method']=='weighted_rrr' and r['bank']=='complete' and float(r['rho'])==.2 and int(r['k'])==int(r['d'])) or (r['bank'] in ('global_trace','global_products') and int(r['k'])==int(r['d']))):continue
    text.append('| '+' | '.join([r['design'],r['amplitude'],r['d'],r['n'],r['bank'],r['method'],r['k'],r['rho'],f"{float(r['mean_loss']):.4f}",f"{float(r['mcse_loss']):.4f}"])+' |')
(dest/'RESULTS.md').write_text('\n'.join(text)+'\n')
print(json.dumps({'figures':10,'reporting':str(dest)}))
