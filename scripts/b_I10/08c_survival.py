"""I10 step 8c: survival (Cox OS ~ footprint + age + late stage + purity per cohort; DL meta). Same code as the survival block of 08_discovery.py, run standalone."""

import sys, numpy as np, pandas as pd, warnings
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import dl_meta
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
S = pd.read_parquet(CA + 'cptac_scores_tc_padj.parquet')
S = S[~S.footprint.str.contains(r'\|')]
F = S.pivot_table(index='sample', columns='footprint', values='nz')
KZ = cl.activity('tc_padj', 'corev3_zcore'); KN = cl.activity('tc_padj', 'corev3_n')
tum = cl.tumour_ids()
kok = ((KN[tum] >= 5).mean(1) >= 0.8); KZ = KZ.loc[kok, tum].T
print('footprints', F.shape[1], 'kinases', KZ.shape[1])
clin = cl.clinical().assign(sample=lambda d: d.case_id + '_T').set_index('sample')
ANN = pd.concat([cl.table(c, 'annotations').assign(sample=lambda d: d.case_id + '_T') for c in cl.COHORTS]).set_index('sample')

# survival
rows = []
for coh in cl.COHORTS:
    ids = cl.tumour_ids(coh); c = clin.reindex(ids)
    T = c.OS_days.astype(float); E = c.OS_event.astype(float)
    if E.sum() < 15: continue
    late = c.Stage.map({'Stage I': 0, 'Stage II': 0, 'Stage III': 1, 'Stage IV': 1}).astype(float)
    pur = c.WES_purity.astype(float); age = c.Age.astype(float)
    for fp in F.columns:
        x = F[fp].reindex(ids); x = (x - x.mean()) / x.std()
        cov = pd.DataFrame({'x': x, 'age': age.fillna(age.median()), 'pur': pur.fillna(pur.median())})
        if late.notna().sum() > 0.5 * len(ids): cov['late'] = late.fillna(0)
        ok = cov.notna().all(1) & T.notna() & E.notna() & (T > 0)
        try:
            r = PHReg(T[ok].values, cov[ok].values, status=E[ok].values).fit()
            rows.append((fp, coh, r.params[0], r.bse[0], int(E[ok].sum()), int(ok.sum())))
        except Exception as e: pass
SV = pd.DataFrame(rows, columns=['footprint', 'cohort', 'logHR', 'se', 'events', 'n'])
SV.to_csv(CA + 'survival_percohort.tsv', sep='\t', index=False)
out = []
for fp, d in SV.groupby('footprint'):
    est, se, p, tau2, I2 = dl_meta(d.logHR.values, d.se.values ** 2)
    out.append(dict(footprint=fp, n_cohorts=len(d), logHR=est, se=se, p=p, I2=I2, n_nominal=int((2 * norm.sf(abs(d.logHR / d.se)) < 0.05).sum()),
                    cohorts=';'.join(f'{c}:{h:+.2f}' for c, h in zip(d.cohort, d.logHR))))
SVm = pd.DataFrame(out); SVm['q'] = multipletests(SVm.p, method='fdr_bh')[1]
SVm.sort_values('p').to_csv(RT + 'T10_survival.tsv', sep='\t', index=False)
print(SVm.sort_values('p').head(15).round(3).to_string())
