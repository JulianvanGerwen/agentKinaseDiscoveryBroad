import sys,numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'/workspace/scripts/infra_compendium'); sys.path.insert(0,'/workspace/scripts/infra_cptac'); import load, cptac_load as cl
T='/workspace/results/val_1/tables/'; F='/workspace/results/val_1/figures/'
R=pd.read_csv(T+'V1_inhibitor_effect.tsv',sep='\t'); Cr=pd.read_csv(T+'V1_activity_corr.tsv',sep='\t')
fig,ax=plt.subplots(1,2,figsize=(11,4.5))
r=R[R.group.isin(['own_inhibitor(primary+profiled)','upstream_MEKi','stimulus_TNF'])&R.effect.notna()].copy(); r['lab']=r.site.str.split('_').str[1].radd(r.kinase+' ')+' ['+r.group.str.replace('own_inhibitor(primary+profiled)','own KI').str.replace('upstream_MEKi','MEKi').str.replace('stimulus_','')+f'] n='+r.n_group.astype(int).astype(str)
ax[0].barh(r.lab,r.effect,color=['#c44' if e<0 else '#48a' for e in r.effect]); ax[0].axvline(0,c='k',lw=.6); ax[0].set_xlabel('median log2FC, group minus all other conditions'); ax[0].tick_params(labelsize=7)
c=Cr[Cr.activity.isin(['loso','core_masked'])&Cr.rho.notna()]; pv=c.pivot(index='site',columns='activity',values='rho')
pv.plot.barh(ax=ax[1],color=['#48a','#c84']); ax[1].set_xlabel('Spearman rho, site vs kinase activity (all conditions)'); ax[1].legend(['core-v3 z (own protein masked)','PSP mean-z leave-self-out'],fontsize=7); ax[1].tick_params(labelsize=7)
plt.tight_layout(); plt.savefig(F+'V1_reporters_compendium.png',dpi=150); plt.close()
M=pd.read_csv(T+'V2_meta.tsv',sep='\t'); M=M[M.k_cohorts>=5]; M=M[M.y.isin(['Hippo_full(tc)','Hippo_noYAP1(masked)'])]
fig,ax=plt.subplots(figsize=(6,3.5)); M=M.reset_index(drop=True); yy=np.arange(len(M))
ax.errorbar(M.meta_r,yy,xerr=[M.meta_r-M.ci_lo,M.ci_hi-M.meta_r],fmt='o',color='#333'); ax.set_yticks(yy); ax.set_yticklabels(M.x+' -> '+M.y,fontsize=7); ax.axvline(0,c='k',lw=.6)
ax.set_xlabel('meta partial Spearman r (purity+proliferation adj), 10 CPTAC cohorts\nhypothesis predicts r<0'); plt.tight_layout(); plt.savefig(F+'V2_hippo_meta.png',dpi=150); plt.close()
C=load.conditions().set_index('conditionID'); X=load.sites_matrix(); y=X.loc['O96013_S148'].dropna(); d=C.loc[y.index].assign(lfc=y)
def g(r):
    a=str(r.agents)
    if a=='PF3758309': return 'PF-3758309'
    if a=='FRAX486': return 'FRAX486'
    if r.study_pmid=='31959955': return 'other Hijazi KI'
    if r.perturbation_class=='kinase_inhibitor': return 'other KI (non-Hijazi)'
    return 'non-inhibitor'
d['grp']=d.apply(g,axis=1); order=['PF-3758309','FRAX486','other Hijazi KI','other KI (non-Hijazi)','non-inhibitor']
fig,ax=plt.subplots(figsize=(6,3.5)); rng=np.random.default_rng(0)
for i,o in enumerate(order):
    v=d.lfc[d.grp==o]; ax.scatter(i+rng.uniform(-.2,.2,len(v)),v,s=8,alpha=.6)
ax.set_xticks(range(5)); ax.set_xticklabels([f'{o}\n(n={(d.grp==o).sum()})' for o in order],fontsize=6); ax.set_ylabel('PAK4 S148 log2FC'); plt.tight_layout(); plt.savefig(F+'V3_pak4_s148.png',dpi=150)
print(d.groupby('grp').lfc.describe().round(2))
n=cl.activity('tc','corev3_n'); print(n.loc[['MET','ALK','LATS1','LATS2','STK38','STK38L']].median(axis=1))
