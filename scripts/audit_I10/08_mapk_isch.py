"""audit_I10 step 8: does the MEK-footprint axis keep its MAPK genetics/kinase meaning after removing the cold-ischaemia component?
Per cohort (tc_padj, own scores): MAPK output = -MEK score; residualised on ischaemia-footprint score (OLS).
Report Spearman with MAPK1/3 core-v3 (tc_padj), RAS/RAF-hotspot AUROC (raw vs residual), and the same for the ischaemia score itself.
Output: A8_mapk_vs_ischaemia.tsv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from common import *
ISCH = ischaemia_set(); MEK = fp_set('KI_MEK'); KA = cl.activity('tc_padj', 'corev3_z')
def hot(coh):
    m = cl.table(coh, 'mutations'); pos = m.protein_change.fillna('').str.extract(r'p\.[A-Z](\d+)')[0].astype(float)
    ras = (m.gene.isin(['KRAS', 'NRAS', 'HRAS']) & pos.isin([12, 13, 61, 117, 146])) | ((m.gene == 'BRAF') & pos.isin([600, 601, 469, 464, 466, 597]))
    return set(m.case_id[ras] + '_T')
rows, pooled = [], []
for coh in COH:
    Z = standardise(cohort(coh))
    D = pd.DataFrame({'mapk': -score(Z, MEK), 'isch': score(Z, ISCH), 'kin': KA.loc[['MAPK1', 'MAPK3']].mean(axis=0).reindex(Z.columns)}).dropna()
    b = np.polyfit(D.isch, D.mapk, 1); D['mapk_res'] = D.mapk - np.polyval(b, D.isch)
    b2 = np.polyfit(D.isch, D.kin, 1); D['kin_res'] = D.kin - np.polyval(b2, D.isch)
    y = D.index.isin(hot(coh)); D['y'] = y
    r = dict(cohort=coh, n=len(D), n_hot=int(y.sum()), rho_mapk_isch=spearmanr(D.mapk, D.isch)[0], rho_kin_isch=spearmanr(D.kin, D.isch)[0],
             rho_mapk_kin=spearmanr(D.mapk, D.kin)[0], rho_mapkres_kinres=spearmanr(D.mapk_res, D.kin_res)[0])
    if 4 <= y.sum() <= len(y) - 4:
        for c in ['mapk', 'mapk_res', 'isch', 'kin', 'kin_res']: r[f'auroc_{c}'] = roc_auc_score(y, D[c])
        pooled.append(D.assign(coh=coh))
    rows.append(r)
P = pd.concat(pooled); out = {'cohort': 'POOLED(within-cohort pct)', 'n_hot': int(P.y.sum())}
for c in ['mapk', 'mapk_res', 'isch', 'kin', 'kin_res']: out[f'auroc_{c}'] = roc_auc_score(P.y, P.groupby('coh')[c].rank(pct=True))
rows.append(out)
T = pd.DataFrame(rows); T.to_csv(RA + 'A8_mapk_vs_ischaemia.tsv', sep='\t', index=False); pd.set_option('display.width', 250); print(T.round(3).to_string())
