import sys, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/scripts/infra_compendium'); import load
from scipy.stats import mannwhitneyu, spearmanr, binomtest
OUT='/workspace/results/val_1/tables/'
C=load.conditions().set_index('conditionID'); S=load.sites().set_index('site')
X=load.sites_matrix()
sign=np.where(C.contrast_reversed.astype(bool),-1.0,1.0); sg=pd.Series(sign,index=C.index)
ok=C.qc_flag.eq('ok')
SITES={ # site: (kinase, expected sign under own inhibitor, extra 'upstream' targets, stimulus)
'P00519_Y185':('ABL1',-1,[],None),'Q14004_T588':('CDK13',-1,[],None),'Q05397_S840':('PTK2',-1,[],None),
'Q05397_Y576':('PTK2',-1,[],None),'Q13177_S197':('PAK2',-1,[],None),
'O14920_S675':('IKBKB',-1,[],'TNF'),'O14920_S697':('IKBKB',-1,[],'TNF'),
'P28482_Y187':('MAPK1',-1,['MAP2K1','MAP2K2'],None),'P27361_Y204':('MAPK3',-1,['MAP2K1','MAP2K2'],None)}
def tg(s): return set(str(s).split(';')) if isinstance(s,str) else set()
C['prim']=C.inhibitor_primary_targets.map(tg); C['prof']=C.inhibitor_profiled_targets.map(tg)
isKI=C.perturbation_class.eq('kinase_inhibitor')&~C.is_combination.astype(bool)
# robust per-condition z for leave-self-out PSP activity
Xv=X.values.astype(float); med=np.nanmedian(Xv,0); mad=np.nanmedian(np.abs(Xv-med),0)*1.4826
Z=pd.DataFrame((Xv-med)/mad,index=X.index,columns=X.columns)
Z=Z*sg.reindex(X.columns).values  # orient
def psp_loso(k, own_uniprot):
    sub=[s for s,v in S.psp_kinases_human.dropna().items() if k in v.split(';') and not s.startswith(own_uniprot+'_')]
    sub=[s for s in sub if s in Z.index]
    M=Z.loc[sub]; n=M.notna().sum(); a=M.mean()*np.sqrt(n); a[n<5]=np.nan
    return a, len(sub)
Acore=load.activity('corev3','z'); Acz=load.activity('corev3','z_core') if False else None
L=pd.read_parquet('/workspace/results/val_1/tables/corev3_masked_long.parquet'); Lc=L; print(L.columns.tolist())
zc=Lc.pivot(index='kinase',columns='condition',values='z_core'); nc=Lc.pivot(index='kinase',columns='condition',values='n_substrates')
Apsp=load.activity('psp_z','z')
rows=[];crow=[]
for site,(k,exp,up,stim) in SITES.items():
    y=(X.loc[site]*sg.reindex(X.columns)).dropna(); y=y[y.index.map(ok)]
    own=pd.Series([bool(({k}&(C.at[c,'prim']|C.at[c,'prof']))) and bool(isKI[c]) for c in y.index],index=y.index)
    own_prim=pd.Series([bool(({k}&C.at[c,'prim'])) and bool(isKI[c]) for c in y.index],index=y.index)
    ups=pd.Series([bool(set(up)&C.at[c,'prim']) and bool(isKI[c]) for c in y.index],index=y.index)
    st=pd.Series([(stim is not None and C.at[c,'stimulus']==stim and C.at[c,'perturbation_class'] in('cytokine','other')) for c in y.index],index=y.index)
    groups={'own_inhibitor(primary+profiled)':own,'own_inhibitor(primary only)':own_prim}
    if up: groups['upstream_MEKi']=ups
    if stim: groups[f'stimulus_{stim}']=st
    for g,m in groups.items():
        rest=~(own|ups|st|own_prim)
        a=y[m]; b=y[rest]
        if len(a)<3: rows.append(dict(site=site,kinase=k,group=g,n_group=len(a),n_other=len(b))); continue
        u=mannwhitneyu(a,b); auc_dec=u.statistic/(len(a)*len(b))  # P(own > other)
        auc = 1-auc_dec if (g.startswith('stimulus')==False) else auc_dec  # AUROC for expected direction: inhibitor expected lower; stimulus expected higher
        # within-study
        st_d=[]
        for s,idx in C.loc[y.index].groupby('study_pmid').groups.items():
            ia=[i for i in idx if m[i]]; ib=[i for i in idx if rest[i]]
            if ia and ib: st_d.append(y[ia].median()-y[ib].median())
        st_d=np.array(st_d); expdir = 1 if g.startswith('stimulus') else -1
        nfav=int((np.sign(st_d)==expdir).sum()) if len(st_d) else 0
        rows.append(dict(site=site,kinase=k,group=g,n_group=len(a),n_studies_group=C.loc[a.index,'study_pmid'].nunique(),n_other=len(b),
          median_group=a.median(),median_other=b.median(),effect=a.median()-b.median(),AUROC_expected_dir=auc,p_MWU=u.pvalue,
          frac_group_in_expdir=float((np.sign(a)==expdir).mean()),n_within_study=len(st_d),within_study_median_delta=np.median(st_d) if len(st_d) else np.nan,
          within_study_n_expdir=nfav,within_study_sign_p=binomtest(nfav,len(st_d),0.5).pvalue if len(st_d) else np.nan))
    # activity correlation (all conditions measured, qc ok, un-flipped orientation is fine since both flipped? use flipped y with Z flipped consistently)
    ap,nsub=psp_loso(k,site.split('_')[0])
    d=pd.DataFrame({'y':y,'loso':ap.reindex(y.index),'core_masked':zc.loc[k].reindex(y.index) if k in zc.index else np.nan,
                    'ncore':nc.loc[k].reindex(y.index) if k in nc.index else np.nan,'psp':Apsp.loc[k].reindex(y.index) if k in Apsp.index else np.nan,
                    'study':C.loc[y.index,'study_pmid']})
    d.loc[d.ncore<1,'core_masked']=np.nan
    for col in ['loso','core_masked','psp']:
        dd=d[['y',col,'study']].dropna()
        if len(dd)<10: crow.append(dict(site=site,kinase=k,activity=col,n=len(dd),n_substrates_psp=nsub)); continue
        r,p=spearmanr(dd.y,dd[col])
        # within-study centred
        c=dd.groupby('study')[['y',col]].transform(lambda v:v-v.mean()); 
        multi=dd.groupby('study').y.transform('size')>1
        rc,pc=spearmanr(c[multi].y,c[multi][col]) if multi.sum()>10 else (np.nan,np.nan)
        crow.append(dict(site=site,kinase=k,activity=col,n=len(dd),n_studies=dd.study.nunique(),rho=r,p=p,rho_within_study=rc,p_within=pc,n_within=int(multi.sum()),n_substrates_psp=nsub))
R=pd.DataFrame(rows); Cr=pd.DataFrame(crow)
R.to_csv(OUT+'V1_inhibitor_effect.tsv',sep='\t',index=False); Cr.to_csv(OUT+'V1_activity_corr.tsv',sep='\t',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(R.round(3).to_string()); print(Cr.round(3).to_string())
