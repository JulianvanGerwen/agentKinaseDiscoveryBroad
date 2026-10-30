"""I10 step 8b: follow-up of candidate mimicry associations with stricter covariates (purity + log TMB + E2F hallmark proliferation)
and comparison with kinase activity (core-v3 z_core, tc_padj) for the matched biology. Output: T9b_followup.tsv"""
import sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
import statsmodels.api as sm
from scipy.stats import norm
sys.path.insert(0,'/workspace/scripts/infra_cptac'); import cptac_load as cl
CA="/workspace/cache/b_I10/"; RT="/workspace/results/b_I10/tables/"
S=pd.read_parquet(CA+'cptac_scores_tc_padj.parquet'); F=S.pivot(index='sample',columns='footprint',values='nz')
KZ=cl.activity('tc_padj','corev3_zcore').T
clin=cl.clinical().assign(sample=lambda d:d.case_id+'_T').set_index('sample')
ANN=pd.concat([cl.table(c,'annotations').assign(sample=lambda d:d.case_id+'_T') for c in cl.COHORTS]).set_index('sample')
def mut(coh,g):
    mb=cl.matrix(coh,'mut_binary'); mb.columns=mb.columns+'_T'
    return mb.loc[g].astype(float) if g in mb.index else None
TESTS=[('mut:TP53',['DNA_damage','KI_DDR','replication_stress','mitosis_cellcycle'],['ATM','ATR','CHEK1','CHEK2','CDK1']),
       ('CIN',['DNA_damage','KI_DDR'],['ATM','ATR','CHEK1','CDK1']),
       ('mut:BRAF',['KI_PI3K_mTOR','PI3Kaxis','KI_MEK'],['MTOR','AKT1','MAPK1']),
       ('immune',['phosphatase_inh','TNF','immune_receptor'],['LCK','ZAP70','SYK','JAK1'])]
rows=[]
for feat,fps,kins in TESTS:
    for sc in fps+['kin:'+k for k in kins]:
        for adj in ['purity','purity+TMB+E2F']:
            per=[]
            for coh in cl.COHORTS:
                ids=cl.tumour_ids(coh); c=clin.reindex(ids)
                if feat.startswith('mut:'):
                    x=mut(coh,feat[4:]); 
                    if x is None: continue
                    x=x.reindex(ids)
                    if x.sum()<5: continue
                elif feat=='CIN': x=c.CNV_index_CIN_index.astype(float)
                else: x=ANN.reindex(ids).ESTIMATE_ImmuneScore.astype(float); x=(x-x.mean())/x.std()
                if feat=='CIN': x=(x-x.mean())/x.std()
                y=KZ[sc[4:]].reindex(ids) if sc.startswith('kin:') else F[sc].reindex(ids)
                D=pd.DataFrame({'x':x,'y':y,'pur':c.WES_purity.astype(float)})
                if adj!='purity': D['tmb']=np.log1p(c.TMB.astype(float)); D['e2f']=ANN.reindex(ids).HALLMARK_E2F_TARGETS.astype(float)
                D=D.fillna(D.median()) if False else D; D=D.apply(lambda s:s.fillna(s.median()) if s.name!='y' and s.name!='x' else s).dropna()
                if len(D)<30: continue
                D['y']=(D.y-D.y.mean())/D.y.std()
                r=sm.OLS(D.y,sm.add_constant(D.drop(columns='y'))).fit()
                per.append((coh,r.params['x'],r.pvalues['x'],len(D)))
            if not per: continue
            P=pd.DataFrame(per,columns=['coh','b','p','n']); z=norm.isf(P.p/2)*np.sign(P.b)
            mz=(z*np.sqrt(P.n)).sum()/np.sqrt(P.n.sum())
            rows.append(dict(feature=feat,score=sc,adjust=adj,meta_z=mz,n_coh=len(P),n_rep=int(((P.p<0.05)&(np.sign(P.b)==np.sign(mz))).sum()),
                             cohorts=';'.join(f'{a}:{b:+.2f}' for a,b in zip(P.coh,P.b))))
T=pd.DataFrame(rows); T.to_csv(RT+'T9b_followup.tsv',sep='\t',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',120); print(T.round(2).to_string(index=False))
