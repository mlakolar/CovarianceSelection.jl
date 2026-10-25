"""Paired contrasts and compact scientific figures; execute under Slurm only."""
import csv,json,math,sys
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);dest=root/'reporting';dest.mkdir(exist_ok=True)
selection=json.loads((root/'grid/selection.json').read_text())
allrows=[];pairrows=[];receipts={}
for stage in ('grid','paired','corrections','confirmation','double_training'):
    folder=root/stage
    if not (folder/'receipt.json').exists():continue
    receipts[stage]=json.loads((folder/'receipt.json').read_text())
    summary=list(csv.DictReader((folder/'summary.csv').open()))
    allrows += [dict(stage=stage,**r) for r in summary]
    if stage=='grid':continue
    sums=defaultdict(lambda:defaultdict(float))
    for path in sorted(folder.glob('part-*.csv')):
        records=defaultdict(dict)
        for r in csv.DictReader(path.open()):
            cid=int(r['config'])%100 if stage=='corrections' else int(r['config'])
            ident=(r['scenario'],r['replicate'],cid,r['point'],r['target'])
            arm=r['adjustment'] if stage=='corrections' else r['index']
            records[ident][arm]=r
        baseline='basic' if stage=='corrections' else 'known'
        for ident,arms in records.items():
            if baseline not in arms:raise ValueError(f'missing baseline {stage} {ident}')
            b=arms[baseline]
            for arm,r in arms.items():
                if arm==baseline:continue
                key=(stage,r['design'],r['n'],r['ntrain'],ident[2],r['point'],r['target'],arm)
                s=sums[key];s['replicates']+=1
                diff=int(r['covered'])-int(b['covered']);s['diff']+=diff;s['diff2']+=diff*diff
                s['baseline_covered']+=int(b['covered']);s['arm_covered']+=int(r['covered'])
                both=r['status']==b['status']=='ok';s['both_successful']+=both
                if both:s['length_diff']+=2*1.959963984540054*(float(r['se'])-float(b['se']))
    for key,s in sorted(sums.items()):
        n=s['replicates'];r=dict(zip(('stage','design','n','ntrain','config','point','target','comparison'),key))
        r.update(replicates=int(n),baseline_coverage=s['baseline_covered']/n,comparison_coverage=s['arm_covered']/n,coverage_difference=s['diff']/n,coverage_difference_mcse=math.sqrt(max(0,(s['diff2']-s['diff']**2/n)/(n-1)/n)),both_successful=int(s['both_successful']),length_difference=s['length_diff']/s['both_successful'] if s['both_successful'] else math.nan)
        pairrows.append(r)
with (dest/'paired.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(pairrows[0]));w.writeheader();w.writerows(pairrows)
with (dest/'summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
# Show independent confirmation only, with frozen development-best specifications.
fig,axes=plt.subplots(3,2,figsize=(10,10),sharex=True)
labels={'opposed_scale':('diagonal_plus','edge'),'noncommuting':('diagonal','varying_edge'),'correlation':('diagonal_1','edge')}
for i,(family,targets) in enumerate(labels.items()):
    cid=selection['designs'][family]['development_best']
    for j,target in enumerate(targets):
        ax=axes[i,j];minimum=1.
        for index,label,color in [('known','Known beta','C0'),('estimated','Estimated beta; independent training n','C1')]:
            rr=sorted([r for r in allrows if r['stage']=='confirmation' and r['design']==family and int(r['config'])==cid and float(r['point'])==0 and r['target']==target and r['index']==index],key=lambda r:int(r['n']))
            if not rr:raise ValueError(f'missing confirmation {family} {target} {index}')
            y=[float(r['coverage']) for r in rr];minimum=min(minimum,min(y))
            ax.errorbar([int(r['n']) for r in rr],y,yerr=[1.96*math.sqrt(v*(1-v)/int(r['replicates'])) for v,r in zip(y,rr)],marker='o',capsize=3,label=label,color=color)
        ax.axhline(.95,color='black',ls='--',lw=1);ax.set_ylim(max(0,min(.85,minimum-.04)),1.005)
        ax.set_title(f'{family}: {target}, t=0');ax.set_xscale('log',base=2);ax.set_xticks([200,400,800,1600],['200','400','800','1600']);ax.grid(alpha=.2)
        if j==0:ax.set_ylabel('Pointwise coverage')
        if i==2:ax.set_xlabel('Precision sample size n (estimated-beta total = 2n)')
handles,names=axes[0,0].get_legend_handles_labels();fig.legend(handles,names,loc='lower center',ncol=2)
fig.suptitle('Basic refitted square-root Lasso; 2000 independent confirmation repetitions')
fig.tight_layout(rect=(0,.04,1,.97))
for ext in ('png','pdf'):fig.savefig(dest/f'central_coverage.{ext}',dpi=160)
plt.close(fig)
# Nonzero points are essential: the correlation family has an exact t=0 symmetry.
for role in ('diagonal','edge'):
    fig,axes=plt.subplots(3,3,figsize=(13,10),sharex=True)
    for i,(family,targets) in enumerate(labels.items()):
        cid=selection['designs'][family]['development_best'];target=targets[0 if role=='diagonal' else 1]
        for j,point in enumerate((-.5,0.,.5)):
            ax=axes[i,j];minimum=1.
            for index,label,color in [('known','Known beta','C0'),('estimated','Estimated beta','C1')]:
                rr=sorted([r for r in allrows if r['stage']=='confirmation' and r['design']==family and int(r['config'])==cid and float(r['point'])==point and r['target']==target and r['index']==index],key=lambda r:int(r['n']))
                y=[float(r['coverage']) for r in rr];minimum=min(minimum,min(y))
                ax.errorbar([int(r['n']) for r in rr],y,yerr=[1.96*math.sqrt(v*(1-v)/int(r['replicates'])) for v,r in zip(y,rr)],marker='o',capsize=3,label=label,color=color)
            ax.axhline(.95,color='black',ls='--',lw=1);ax.set_ylim(max(0,min(.85,minimum-.04)),1.005)
            ax.set_title(f'{family}, t={point:g}');ax.set_xscale('log',base=2);ax.set_xticks([200,400,800,1600],['200','400','800','1600']);ax.grid(alpha=.2)
            if j==0:ax.set_ylabel('Pointwise coverage')
            if i==2:ax.set_xlabel('Precision n; estimated-beta total = 2n')
    handles,names=axes[0,0].get_legend_handles_labels();fig.legend(handles,names,loc='lower center',ncol=2)
    fig.suptitle(f'{role.capitalize()} coverage; frozen basic procedure; 2000 repetitions')
    fig.tight_layout(rect=(0,.04,1,.97))
    for ext in ('png','pdf'):fig.savefig(dest/f'{role}_coverage.{ext}',dpi=160)
    plt.close(fig)
receipt={'stages':receipts,'paired_rows':len(pairrows),'summary_rows':len(allrows),'confirmation_replicates':2000,'primary_corrections':False,'full_failures_retained':True,'figures':['central_coverage.png','central_coverage.pdf','diagonal_coverage.png','diagonal_coverage.pdf','edge_coverage.png','edge_coverage.pdf']}
(dest/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
# Compact values for interpretation. All numerical aggregation occurs here under Slurm.
for family,targets in labels.items():
    cid=selection['designs'][family]['development_best']
    print('\n',family,'frozen_config',cid)
    for n in (200,400,800,1600):
        cells=[r for r in allrows if r['stage']=='confirmation' and r['design']==family and int(r['config'])==cid and float(r['point'])==0 and int(r['n'])==n and r['target'] in targets]
        print(n,[(r['target'],r['index'],round(float(r['coverage']),4),round(float(r['bias']),4),round(float(r['empirical_sd']),4),round(float(r['rms_se']),4)) for r in cells])
for family,targets in labels.items():
    cid=selection['designs'][family]['development_best']
    print('ALL_POINTS',family,'config',cid)
    for n in (200,400,800,1600):
        print(n,[(r['point'],r['target'],r['index'],round(float(r['coverage']),4)) for r in allrows if r['stage']=='confirmation' and r['design']==family and int(r['config'])==cid and int(r['n'])==n and r['target'] in targets])
print(json.dumps(receipt))
