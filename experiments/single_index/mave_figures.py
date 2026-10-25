"""Scientific index-estimation figure from validated Slurm summaries."""
import csv,sys,json,os
for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
 os.environ[name]="1"
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);out=root/'interpretation';out.mkdir(exist_ok=True)
main=list(csv.DictReader((root/'confirmation/summary.csv').open()))
hybrid=list(csv.DictReader((root/'hybrid/summary.csv').open()))
series=[('stein_screen10','Stein screened/refit',main,'#1f77b4','o',None),
        ('one_update','One dense MAVE update',main,'#17becf','s',None),
        ('full_mave','Iterated dense MAVE (cap 200)',main,'#d62728','^',None),
        ('mave_select_rrr_refit','MAVE selection + RRR refit',hybrid,'#9467bd','D',None),
        ('one_adaptive_select_rrr_refit','One sparse update + RRR refit',hybrid,'#ff7f0e','v',None),
        ('stein_threshold_rrr_refit','Stein threshold + RRR refit (C=1)',hybrid,'#2ca02c','P',1.),
        ('stein_threshold_rrr_refit','Stein threshold + RRR refit (C=1.5)',hybrid,'#8c564b','*',1.5),
        ('oracle_support_rrr','True-support RRR (infeasible)',hybrid,'#333333','x',None)]
names={'correlation':'Correlation','opposed_scale':'Opposed scale','fixed_diagonal':'Fixed precision diagonal','noncommuting':'Noncommuting'}
fig,axes=plt.subplots(2,2,figsize=(11.5,8),sharex=True,sharey=True)
for ax,(design,title) in zip(axes.flat,names.items()):
 for method,label,rows,color,marker,level in series:
  selected=[r for r in rows if r['design']==design and int(r['p'])==40 and r['method']==method and (level is None or float(r['level'])==level)]
  selected.sort(key=lambda r:int(r['n']));assert len(selected)==4,(design,method)
  x=[int(r['n']) for r in selected];y=[float(r['mean_loss']) for r in selected];se=[float(r['mcse_loss']) for r in selected]
  ax.errorbar(x,y,yerr=se,label=label,color=color,marker=marker,lw=1.3,ms=4,capsize=2,ls='--' if method=='oracle_support_rrr' else '-')
 ax.set_title(title,fontsize=11);ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks([200,400,800,1600],['200','400','800','1600']);ax.set_ylim(.0003,1.3);ax.grid(True,alpha=.2)
for ax in axes[:,0]:ax.set_ylabel('Mean squared-sine direction loss')
for ax in axes[-1,:]:ax.set_xlabel('Index-training observations')
handles,labels=axes.flat[0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',ncol=2,fontsize=8,bbox_to_anchor=(.5,.008))
fig.suptitle('Index estimation: 20 Gaussian covariates, 40 responses',fontsize=13)
fig.text(.5,.155,'300 paired replications per setting; bars show one Monte Carlo standard error.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.17,1,.95))
fig.savefig(out/'index-comparison.png',dpi=180);fig.savefig(out/'index-comparison.pdf');plt.close(fig)

(out/'figure-environment.json').write_text(json.dumps(dict(python=sys.version,matplotlib=matplotlib.__version__,source='Validated confirmation and hybrid summary tables; no refitting.'),indent=2)+'\n')
