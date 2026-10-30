import os,sys,numpy as np,pandas as pd
os.environ['OMP_NUM_THREADS']='4'
sys.path.insert(0,'/workspace/scripts/infra_cptac'); sys.path.insert(0,'/workspace/scripts/infra_compendium/kinact_mvp_pkg')
import cptac_load as cl
from kinact_mvp import KinaseActivityModel
from threadpoolctl import threadpool_limits
from scipy.stats import rankdata, t as tdist, norm, mannwhitneyu
OUT='/workspace/results/val_1/tables/'
model=KinaseActivityModel.from_config('/workspace/results/infra_cptac/model_corev3/config.json')
MASK=['P46937','P08581','Q9UM73','P28482','P27361','O95835','Q9NRM7','Q15208','Q9Y2H1']
HIP=['LATS1','LATS2','STK38','STK38L']
S=cl.samples().set_index('sample_id')
def rz(df): # row-wise robust? use per-kinase rank-normal not needed; z-score each kinase across tumours
    return (df.sub(df.mean(1),axis=0)).div(df.std(1),axis=0)
def partial_spearman(x,y,cov):
    d=pd.concat([x,y,cov],axis=1).dropna()
    if len(d)<30: return np.nan,len(d)
    R=d.apply(rankdata).values; R=(R-R.mean(0))/R.std(0)
    C=np.c_[np.ones(len(R)),R[:,2:]]
    rx=R[:,0]-C@np.linalg.lstsq(C,R[:,0],rcond=None)[0]; ry=R[:,1]-C@np.linalg.lstsq(C,R[:,1],rcond=None)[0]
    return np.corrcoef(rx,ry)[0,1],len(d)
def meta(r,n):
    r=np.array(r);n=np.array(n); m=np.isfinite(r); r,n=r[m],n[m]
    z=np.arctanh(r); v=1/(n-3-3); w=1/v; zf=(w*z).sum()/w.sum(); Q=(w*(z-zf)**2).sum(); k=len(z)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); ws=1/(v+tau2); zr=(ws*z).sum()/ws.sum(); se=np.sqrt(1/ws.sum())
    p=2*norm.sf(abs(zr/se)); return np.tanh(zr),np.tanh(zr-1.96*se),np.tanh(zr+1.96*se),p,k,int((r<0).sum())
rows=[];grp=[]
for c in cl.COHORTS:
    t=cl.tumour_ids(c); V=cl.matrix(c,'input_tc')[t]
    Vm=V.loc[~V.index.str.split('_').str[0].isin(MASK)]
    with threadpool_limits(limits=4):
        rf=model.score(V,n_jobs=4); rm=model.score(Vm,n_jobs=4)
    zf=rf.z_core.T if hasattr(rf,'z_core') else None
    Af=rf.z.T; Am=rm.z.T; Nm=rm.n_substrates.T
    ann=cl.table(c,'annotations').set_index('case_id'); clin=cl.table(c,'clinical').set_index('case_id')
    cases=S.loc[t,'case_id']
    pur=clin['WES_purity'].reindex(cases.values).values; 
    est=ann['ESTIMATE_ESTIMATEScore'].reindex(cases.values).values
    pur=pd.Series(pur,index=t); est=pd.Series(est,index=t)
    pur=pur.fillna(-est.rank(pct=True)*0) if False else pur
    purity=pur.where(pur.notna(), np.nan)
    if purity.isna().mean()>0.3: purity=-est   # fallback
    prolif=pd.Series(((ann['HALLMARK_G2M_CHECKPOINT']+ann['HALLMARK_E2F_TARGETS'])/2).reindex(cases.values).values,index=t)
    cov=pd.concat([purity.rename('pur'),prolif.rename('prolif')],axis=1)
    pad=cl.matrix(c,'input_tc_padj'); ph=cl.matrix(c,'input_tc')
    yap=pad.loc['P46937_S127',t] if 'P46937_S127' in pad.index else pd.Series(np.nan,index=t)
    ax=lambda A:{'MET':A.loc['MET'],'ALK':A.loc['ALK'],'RTK_MAPK':rz(A.loc[['EGFR','MET','ALK','ERBB2','MAPK1','MAPK3']].dropna(how='all')).mean(0),'MAPK1_3':rz(A.loc[['MAPK1','MAPK3']]).mean(0)}
    Xs={'full':ax(Af),'clean':ax(Am)}
    hip_full=rz(Af.loc[HIP]).mean(0); hip_clean=rz(Am.loc[HIP]).mean(0)
    Ys={'Hippo_full(tc)':hip_full,'Hippo_noYAP1(masked)':hip_clean,'YAP1_S127(padj)':yap}
    for k in HIP: Ys[k+'_noYAP1']=Am.loc[k]
    for xn in ['MET','ALK','RTK_MAPK','MAPK1_3']:
        for yn,y in Ys.items():
            # predictors from clean version for masked outcome, full otherwise
            x=Xs['clean'][xn] if 'noYAP1' in yn or yn.startswith('YAP1') else Xs['full'][xn]
            xx=x if not (xn in('MET','ALK') and (Nm.loc[xn]<1).all()) else x
            r,n=partial_spearman(x.reindex(t),y.reindex(t),cov); r0,_=partial_spearman(x.reindex(t),y.reindex(t),pd.DataFrame(index=t,columns=[]).assign(c=0.0)) if False else (np.nan,0)
            rows.append(dict(cohort=c,x=xn,y=yn,r_partial=r,n=n))
    # MET-high
    metp=cl.matrix(c,'protein'); mp=metp.loc[[i for i in metp.index if i.startswith('P08581')][0]] if any(i.startswith('P08581') for i in metp.index) else None
    try:
        gi=cl.matrix(c,'cnv_gistic'); gm=[i for i in gi.index if i.split('_')[0] in('MET',)] 
    except Exception as e: gi=None
    cid=lambda: None
    if mp is not None:
        mp=mp.reindex(t) if set(t)&set(mp.index) else mp.reindex(t)
    grp.append((c,t,mp,hip_full,hip_clean,yap,cov))
R=pd.DataFrame(rows); R.to_csv(OUT+'V2_partial_corr_by_cohort.tsv',sep='\t',index=False)
M=[]
for (xn,yn),g in R.groupby(['x','y']):
    r,lo,hi,p,k,nneg=meta(g.r_partial,g.n); M.append(dict(x=xn,y=yn,meta_r=r,ci_lo=lo,ci_hi=hi,p=p,k_cohorts=k,k_negative=nneg))
M=pd.DataFrame(M); M.to_csv(OUT+'V2_meta.tsv',sep='\t',index=False); pd.set_option('display.width',200); print(M.round(3).to_string())
import pickle; pickle.dump(grp,open('<SCRATCH>/grp.pkl','wb'))
