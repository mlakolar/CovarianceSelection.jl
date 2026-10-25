"""Paired confirmation comparisons and standalone scientific figures (Slurm only)."""
import csv,json,math,statistics,sys
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(path):
    with path.open() as f:return list(csv.DictReader(f))

def write(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader();w.writerows(rows)

def label(k):
    if k=='validation':return 'Held-out selection'
    if k=='oracle':return 'Oracle among nine'
    bank,L,method,support=k.split('-')
    suffix=bank if L=='0' else f'{bank}, L={L}'
    return f"{method}{' k='+support if method=='screen' else ''}: {suffix}"

def main(root):
    root=Path(root); conf=root/'confirmation_final'; dev=root/'development'
    shortlist=json.loads((root/'experiments/single_index/multiprobe_confirmation_final.shortlist.json').read_text())
    default=shortlist['fixed_default']; ids=[c['id'] for c in shortlist['candidates']]
    rr=read(conf/'summary.csv')
    selected=[]
    for match in (lambda k:'-ols-' in k,lambda k:'-screen-' in k,
                  lambda k:not ('-ols-' in k or '-screen-' in k)):
        selected.append(next(k for k in ids if match(k)))
    selected+=['validation','oracle']
    selected[3:3]=['complete-0-ols-20','complete-0-screen-6']
    figdir=root/'figures';figdir.mkdir(exist_ok=True)
    for metric,ylabel,filename in [('mean_loss','Mean projective error','direction_error'),
                                 ('prob_loss_below_01','Probability(error < 0.1)','reliable_direction')]:
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),sharey=True)
        for ax,design in zip(axes,('chain','pairs')):
            for k in selected:
                rows=sorted([r for r in rr if r['design']==design and r['p']=='40' and
                             r['signal']=='1.0' and r['beta_kind']=='equal' and r['id']==k],
                            key=lambda r:int(r['n']))
                x=[int(r['n']) for r in rows];y=[float(r[metric]) for r in rows]
                if metric=='mean_loss':
                    ax.errorbar(x,y,yerr=[1.96*float(r['se_loss']) for r in rows],
                                marker='o',markersize=3,linewidth=1.3,label=label(k),capsize=2)
                else:ax.plot(x,y,'-o',markersize=3,linewidth=1.3,label=label(k))
            ax.set_title('Chain precision' if design=='chain' else 'Constant-variance pairs')
            ax.set_xscale('log',base=2);ax.set_xticks([100,200,400,800,1600,3200])
            ax.set_xticklabels(['100','200','400','800','1600','3200'],rotation=30)
            ax.set_xlabel('Sample size');ax.set_ylim(-.02,1.02);ax.grid(alpha=.2)
        axes[0].set_ylabel(ylabel)
        handles,labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
        fig.subplots_adjust(bottom=.30,wspace=.08,top=.9)
        fig.savefig(figdir/f'{filename}.pdf',bbox_inches='tight')
        fig.savefig(figdir/f'{filename}.png',dpi=180,bbox_inches='tight');plt.close(fig)

    diagnostic=read(dev/'summary.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True)
    for ax,design in zip(axes,('chain','pairs')):
        for k in ['mixed-60-ols-20','mixed-60-second-20','mixed-60-second_split-20']:
            rows=sorted([r for r in diagnostic if r['even']=='True' and r['design']==design
                         and r['id']==k],key=lambda r:int(r['n']))
            ax.errorbar([int(r['n']) for r in rows],[float(r['mean_loss']) for r in rows],
                        yerr=[1.96*float(r['se_loss']) for r in rows],label=label(k),marker='o',capsize=2)
        ax.set_title(design);ax.set_xlabel('Sample size');ax.set_ylim(0,1);ax.grid(alpha=.2)
    axes[0].set_ylabel('Mean projective error')
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',fontsize=8,ncol=1)
    fig.subplots_adjust(bottom=.3,wspace=.08)
    fig.savefig(figdir/'even_link.pdf',bbox_inches='tight')
    fig.savefig(figdir/'even_link.png',dpi=180,bbox_inches='tight');plt.close(fig)

    # Preserve paired differences rather than subtracting independent SEs.
    data=defaultdict(dict)
    for path in sorted(conf.glob('part-*.csv')):
        for r in read(path):
            key=(r['scenario'],r['replicate'])
            name=r['selected_id'] if r['selector']=='fixed' else r['selector']
            data[key][name]=float(r['loss']) if r['status']=='ok' else 1.0
    differences=defaultdict(list)
    for (sid,rep),values in data.items():
        for k in ids+['validation','oracle']:
            if k!=default:differences[(sid,k)].append(values[k]-values[default])
    paired=[]
    for (sid,k),v in sorted(differences.items()):
        paired.append(dict(scenario=sid,id=k,reference=default,replicates=len(v),
                           mean_difference=statistics.mean(v),
                           se_difference=statistics.stdev(v)/math.sqrt(len(v)),
                           probability_better=sum(x<0 for x in v)/len(v)))
    write(conf/'paired_differences.csv',paired)
    # Exact identification diagnostic for RANDOM PRODUCT probes in the paired
    # design. This uses the generator only to explain results, never for fitting.
    diagnostic=[]
    for p in (40,160):
        total=p*(p-1)//2; informative=p//2
        for L in (10,30,60):
            diagnostic.append(dict(p=p,total_products=total,informative_products=informative,L=L,
                probability_zero_first_order_signal=math.comb(total-informative,L)/math.comb(total,L)))
    write(conf/'random_product_signal.csv',diagnostic)
    selected_rows=[r for r in rr if r['id'] in selected]
    write(conf/'compact_comparison.csv',selected_rows)
    print('Figures and paired differences complete:',default,len(paired))

if __name__=='__main__':main(sys.argv[1])
