"""audit_I10 step 9: correlates of footprint scores (own, tc_padj) with sample-level covariates, per cohort (Spearman; |rho| for PCs).
Output: A9_correlates.tsv (footprint x covariate, median over cohorts + range)"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from common import *
ISCH = ischaemia_set(); FPS = ['KI_MEK', 'EGF', 'GF_all', 'TNF', 'DNA_damage', 'phosphatase_inh', 'GPCR_ligand', 'metformin']
rows = []
for coh in COH:
    Z = standardise(cohort(coh))
    ann = cl.table(coh, 'annotations').assign(s=lambda q: q.case_id + '_T').set_index('s').reindex(Z.columns)
    cli = cl.table(coh, 'clinical').assign(s=lambda q: q.case_id + '_T').set_index('s').reindex(Z.columns)
    A = Z[Z.notna().mean(1) >= 0.5].fillna(0); U, s, Vt = np.linalg.svd(A.values - A.values.mean(1, keepdims=True), full_matrices=False)
    COV = {'purity': cli.WES_purity, 'ESTIMATE_stromal': ann.ESTIMATE_StromalScore, 'ESTIMATE_immune': ann.ESTIMATE_ImmuneScore,
           'E2F_hallmark': ann.HALLMARK_E2F_TARGETS, 'G2M_hallmark': ann.HALLMARK_G2M_CHECKPOINT, 'hypoxia_hallmark': ann.get('HALLMARK_HYPOXIA'),
           'cold_ischaemia_fp_score': score(Z, ISCH), '|PC1|': pd.Series(Vt[0], index=Z.columns), '|PC2|': pd.Series(Vt[1], index=Z.columns),
           '|PC3|': pd.Series(Vt[2], index=Z.columns)}
    for g in FPS:
        sc = score(Z, fp_set(g))
        for k, v in COV.items():
            if v is None: continue
            v = pd.to_numeric(v, errors='coerce'); ok = v.notna() & sc.notna()
            r = spearmanr(sc[ok], v[ok])[0]
            rows.append(dict(cohort=coh, footprint=g, covariate=k, rho=abs(r) if k.startswith('|') else r))
T = pd.DataFrame(rows)
S = T.groupby(['footprint', 'covariate']).rho.agg(['median', 'min', 'max']).round(2)
S['cell'] = S.apply(lambda r: f"{r['median']:+.2f} ({r['min']:+.2f}..{r['max']:+.2f})", axis=1)
W = S.cell.unstack(); W.to_csv(RA + 'A9_correlates.tsv', sep='\t'); T.to_csv(CA + 'A9_correlates_percohort.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); print(S['median'].unstack().to_string())
