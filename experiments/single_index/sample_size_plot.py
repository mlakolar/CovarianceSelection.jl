"""Standalone coverage figures; execute in the existing Discovery report env."""
import csv,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS=['#175d91','#ca6b20','#666666']

def read(path):return list(csv.DictReader(Path(path).open()))

def panel(ax,rows,modes,field,entry):
    for color,(mode,label) in zip(COLORS,modes):
        rr=[r for r in rows if r[field]==mode and (int(r['entry_a']),int(r['entry_b']))==entry]
        nn=sorted({int(r['scenario'].split('_n')[1].split('_')[0]) for r in rr})
        central=[];low=[];high=[];loerr=[];hierr=[]
        for n in nn:
            at=[r for r in rr if int(r['scenario'].split('_n')[1].split('_')[0])==n]
            mid=next(r for r in at if float(r['point'])==0)
            central.append(float(mid['coverage']))
            low.append(min(float(r['coverage']) for r in at));high.append(max(float(r['coverage']) for r in at))
            # Wilson endpoints can differ from 0 or 1 by machine roundoff.
            loerr.append(max(0.,float(mid['coverage'])-float(mid['coverage_lower'])))
            hierr.append(max(0.,float(mid['coverage_upper'])-float(mid['coverage'])))
        ax.plot(nn,central,'o-',color=color,label=label,markersize=3,lw=1.2)
        ax.fill_between(nn,low,high,color=color,alpha=.10)
        ax.errorbar(nn,central,yerr=[loerr,hierr],fmt='none',ecolor=color,lw=.7,capsize=2)
    ax.axhline(.95,color='black',ls='--',lw=.8)
    ax.set_xscale('log',base=2);ax.set_xticks(nn,labels=[str(n) for n in nn],rotation=40)
    ax.set_ylim(0,1.01);ax.grid(alpha=.2);ax.set_xlabel('Sample size')

def plots(root):
    root=Path(root)
    if (root/'known/summary.csv').exists():
        rows=read(root/'known/summary.csv')
        fig,axes=plt.subplots(3,4,figsize=(13,8))
        for j,(prefix,title,entries) in enumerate([
            ('chain_', 'Chain, p=40',[(1,1),(1,2),(1,4)]),
            ('chain_', 'Chain, p=160',[(1,1),(1,2),(1,4)]),
            ('chain_', 'Chain, p=400',[(1,1),(1,2),(1,4)]),
            ('hub6_', 'Hub degree 6, p=160',[(1,1),(1,2),(2,3)])]):
            p=(40,160,400,160)[j]
            rr=[r for r in rows if r['scenario'].startswith(prefix) and r['scenario'].endswith(f'_p{p}')]
            for i,entry in enumerate(entries):
                panel(axes[i,j],rr,[('union_refit','Union refit'),('corrected','Corrections'),('oracle_nuisance','True nuisance')],'variant',entry)
                axes[i,j].set_ylim((.5 if j==3 else .8) if i==0 else .9,1.005)
                if i==0:axes[i,j].set_title(title)
                if j==0:axes[i,j].set_ylabel(('Diagonal','Active edge','Zero edge')[i]+' coverage')
        axes[0,0].legend(fontsize=8,loc='lower right')
        fig.suptitle('Known index: 1,000 replications per setting\nLines: t=0; ribbons: range across t=-0.5,0,0.5; bars: 95% Monte Carlo intervals',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.94));fig.savefig(root/'known/coverage.png',dpi=180);fig.savefig(root/'known/coverage.pdf');plt.close(fig)
    if (root/'estimated/summary.csv').exists():
        rows=read(root/'estimated/summary.csv');fig,axes=plt.subplots(1,3,figsize=(11,3.7),sharey=True)
        for ax,entry,title in zip(axes,[(1,1),(1,2),(1,4)],['Diagonal','Active edge','Zero edge']):
            # Normalize names for the same n parser as the known-index figure.
            rr=[dict(r,scenario=r['scenario'].replace('index_n','index__n')) for r in rows]
            panel(ax,rr,[('known','Known index'),('mean_ols','Selected support + OLS'),('oracle_support_ols','True support + OLS (infeasible)')],'index_method',entry)
            ax.set_title(title)
        axes[0].set_ylabel('Coverage');axes[0].legend(fontsize=8)
        fig.suptitle('Estimated index, p=40: failures count as noncoverage\nLines: t=0; ribbons: range across three locations; bars: Monte Carlo intervals',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.87));fig.savefig(root/'estimated/coverage.png',dpi=180);fig.savefig(root/'estimated/coverage.pdf');plt.close(fig)
    if (root/'bootstrap/summary.csv').exists():
        rows=read(root/'bootstrap/summary.csv');fig,axes=plt.subplots(2,3,figsize=(11,6.5))
        for i,method in enumerate(('mean_ols','oracle_support_ols')):
            rr=[r for r in rows if r['index_method']==method]
            for j,entry in enumerate(((1,1),(1,2),(1,4))):
                panel(axes[i,j],rr,[('wald','Wald'),('fixed_index','Fixed-index bootstrap-t'),('full_pipeline','Full two-stage bootstrap-t')],'interval_method',entry)
                if i==0:axes[i,j].set_title(('Diagonal','Active edge','Zero edge')[j])
                if i==1:axes[i,j].set_ylim(.8,1.005)
            axes[i,0].set_ylabel(('Selected support + OLS','True support + OLS (infeasible)')[i]+'\nCoverage')
        axes[0,0].legend(fontsize=7)
        fig.suptitle('Pairs bootstrap: 300 outer replications, 399 draws\nUnavailable intervals count as noncoverage; lines: t=0; ribbons: three-location range',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.90));fig.savefig(root/'bootstrap/coverage.png',dpi=180);fig.savefig(root/'bootstrap/coverage.pdf');plt.close(fig)

if __name__=='__main__':plots(sys.argv[1])
