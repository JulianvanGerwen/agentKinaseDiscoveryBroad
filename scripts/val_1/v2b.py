import sys,pickle,numpy as np,pandas as pd,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/scripts/infra_cptac'); import cptac_load as cl
from scipy.stats import mannwhitneyu, spearmanr, norm
from scipy.stats import rankdata
src=open('/workspace/scripts/val_1/v2.py').read(); exec(src[src.index('def rz'):src.index('rows=[];grp=[]')])
OUT='/workspace/results/val_1/tables/'
grp=pickle.load(open('<SCRATCH>/grp.pkl','rb'))
rows=[];cov_rows=[];mrows=[]
for c,t,mp,hf,hc,yap,cov in grp:
    pr=cl.matrix(c,'protein'); pi=[i for i in pr.index if i.split('_')[0]=='P08581' or i=='MET']
    ph=cl.matrix(c,'input_tc'); pad=cl.matrix(c,'input_tc_padj'); rawph=cl.matrix(c,'phospho')
    cov_rows.append(dict(cohort=c,n_tumours=len(t),YAP1_S127_tc=int(ph.loc['P46937_S127',t].notna().sum()) if 'P46937_S127' in ph.index else 0,
        YAP1_S127_padj=int(pad.loc['P46937_S127',t].notna().sum()) if 'P46937_S127' in pad.index else 0,
        YAP1_S127_raw=int(rawph.loc['P46937_S127',[x for x in rawph.columns if x in t]].notna().sum()) if 'P46937_S127' in rawph.index else 0,MET_protein=len(pi)))
    ycl=ph.loc['P46937_S127',t] if 'P46937_S127' in ph.index else pd.Series(np.nan,index=t)
    if not pi: continue
    m=pr.loc[pi[0],[x for x in pr.columns if x in t]].reindex(t)
    n=m.notna().sum(); 
    if n<20: continue
    hi=m>=m.quantile(.9); lo=m.notna()&~hi
    # also gistic amp
    for yn,y in [('Hippo_full',hf),('Hippo_noYAP1',hc),('YAP1_S127_tc',ycl)]:
        a=y[hi].dropna(); b=y[lo].dropna()
        if len(a)>=3 and len(b)>=10:
            d=(a.mean()-b.mean())/np.sqrt(((len(a)-1)*a.var()+(len(b)-1)*b.var())/(len(a)+len(b)-2))
            mrows.append(dict(cohort=c,y=yn,n_hi=len(a),n_other=len(b),cohen_d=d,p=mannwhitneyu(a,b).pvalue,var_d=1/len(a)+1/len(b)+d*d/(2*(len(a)+len(b)))))
    # partial r MET protein vs outcomes
    for yn,y in [('Hippo_full',hf),('Hippo_noYAP1',hc),('YAP1_S127_tc',ycl),('YAP1_S127_padj',yap)]:
        r,nn=partial_spearman(m,y.reindex(t),cov); rows.append(dict(cohort=c,x='MET_protein',y=yn,r_partial=r,n=nn))
    # unadjusted RTK/MAPK vs Hippo
    for yn,y in [('Hippo_full',hf)]:
        pass
M=pd.DataFrame(mrows); P=pd.DataFrame(rows); pd.DataFrame(cov_rows).to_csv(OUT+'V2_yap1_coverage.tsv',sep='\t',index=False)
print(pd.DataFrame(cov_rows).to_string())
out=[]
for yn,g in M.groupby('y'):
    w=1/g.var_d; d=(w*g.cohen_d).sum()/w.sum(); Q=(w*(g.cohen_d-d)**2).sum(); k=len(g)
    tau=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); ws=1/(g.var_d+tau); dr=(ws*g.cohen_d).sum()/ws.sum(); se=np.sqrt(1/ws.sum())
    out.append(dict(y=yn,meta_d_METhigh_vs_other=dr,lo=dr-1.96*se,hi=dr+1.96*se,p=2*norm.sf(abs(dr/se)),k=k,k_neg=int((g.cohen_d<0).sum())))
print(pd.DataFrame(out).round(3).to_string()); pd.DataFrame(out).to_csv(OUT+'V2_METhigh_meta.tsv',sep='\t',index=False); M.to_csv(OUT+'V2_METhigh_by_cohort.tsv',sep='\t',index=False)
out=[]
for yn,g in P.groupby('y'):
    r,lo,hi,p,k,kn=meta(g.r_partial,g.n); out.append(dict(x='MET_protein',y=yn,meta_r=r,lo=lo,hi=hi,p=p,k=k,k_neg=kn))
print(pd.DataFrame(out).round(3).to_string()); pd.DataFrame(out).to_csv(OUT+'V2_METprotein_meta.tsv',sep='\t',index=False)
