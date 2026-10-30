"""audit_I10 step 7: is 'immune -> pervanadate mimicry' specific to the PTP-inhibitor footprint?
Per cohort (tc_padj, own scores): OLS z(pervanadate score) ~ z(immune) + purity [+ MEK score + TNF score + EGF score + all-site PC1..3].
Also: same model for the pS/T-only, immune-protein-depleted pervanadate set; site overlap with TNF / EGF / GF_all footprints;
tumour correlation of the pervanadate score with TNF / EGF / MEK scores; pY-site measurability in CPTAC.
Output: A7_pervanadate_specificity.tsv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import statsmodels.api as sm
from scipy.stats import norm
from common import *
SA = site_annot(); PV = fp_set('phosphatase_inh'); HPA = set(open(CA + 'immune_uniprot.txt').read().split())
res = SA.residue.reindex(PV.index).fillna('S').values
PVst = PV[(res != 'Y') & ~np.isin([s.split('_')[0] for s in PV.index], list(HPA))]
ISCH = ischaemia_set()
O = {g: fp_set(g) for g in ['TNF', 'EGF', 'GF_all', 'KI_MEK', 'immune_receptor', 'GPCR_ligand']}
rows = [dict(cohort='ALL', model=f'site overlap with {g}', value=len(PV.index.intersection(s.index))) for g, s in O.items()]
per = {}
for coh in COH:
    Z = standardise(cohort(coh))
    ann = cl.table(coh, 'annotations').assign(s=lambda q: q.case_id + '_T').set_index('s').reindex(Z.columns)
    cli = cl.table(coh, 'clinical').assign(s=lambda q: q.case_id + '_T').set_index('s').reindex(Z.columns)
    pur = cli.WES_purity.astype(float); pur = pur.fillna(pur.median())
    A = Z[Z.notna().mean(1) >= 0.5].fillna(0); U, s, Vt = np.linalg.svd(A.values - A.values.mean(1, keepdims=True), full_matrices=False)
    D = pd.DataFrame({'pv': score(Z, PV), 'pvst': score(Z, PVst), 'imm': ann.ESTIMATE_ImmuneScore.astype(float), 'pur': pur,
                      'mek': score(Z, O['KI_MEK']), 'tnf': score(Z, O['TNF']), 'egf': score(Z, O['EGF']),
                      'isch': score(Z, ISCH), 'pc1': Vt[0], 'pc2': Vt[1], 'pc3': Vt[2]}, index=Z.columns).dropna()
    D = (D - D.mean()) / D.std()
    rows.append(dict(cohort=coh, model='n measured pY sites (>=25% tumours) in cohort', value=int((SA.residue.reindex(Z.index) == 'Y').sum())))
    rows.append(dict(cohort=coh, model='n pervanadate pY sites measured', value=int(((SA.residue.reindex(PV.index) == 'Y') & PV.index.isin(Z.index)).sum())))
    for g in ['tnf', 'egf', 'mek', 'isch']: rows.append(dict(cohort=coh, model=f'r(pervanadate, {g})', value=D.pv.corr(D[g])))
    for y in ['pv', 'pvst']:
        for nm, cov in [('imm+pur', []), ('+MEK+TNF+EGF', ['mek', 'tnf', 'egf']), ('+ischaemia', ['isch']), ('+PC1-3', ['pc1', 'pc2', 'pc3']), ('+all', ['mek', 'tnf', 'egf', 'isch', 'pc1', 'pc2', 'pc3'])]:
            r = sm.OLS(D[y], sm.add_constant(D[['imm', 'pur'] + cov])).fit()
            per.setdefault((y, nm), []).append((r.params['imm'], r.pvalues['imm'], len(D)))
            rows.append(dict(cohort=coh, model=f'{y} ~ immune {nm}', value=r.params['imm'], p=r.pvalues['imm']))
for (y, nm), v in per.items():
    P = pd.DataFrame(v, columns=['b', 'p', 'n']); z = norm.isf(P.p / 2) * np.sign(P.b)
    rows.append(dict(cohort='META', model=f'{y} ~ immune {nm}', value=(z * np.sqrt(P.n)).sum() / np.sqrt(P.n.sum()), p=int(((P.p < .05) & (P.b > 0)).sum())))
T = pd.DataFrame(rows); T.to_csv(RA + 'A7_pervanadate_specificity.tsv', sep='\t', index=False)
pd.set_option('display.width', 200); print(T[T.cohort.isin(['ALL', 'META'])].round(3).to_string())
t = T[~T.cohort.isin(['ALL', 'META'])]; print(t[t.model.str.startswith(('r(', 'n '))].pivot(index='cohort', columns='model', values='value').round(2).to_string())
