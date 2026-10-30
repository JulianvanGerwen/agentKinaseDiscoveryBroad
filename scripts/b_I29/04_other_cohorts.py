"""I29 H6: EGFR-altered vs WT in other CPTAC cohorts (LSCC, HNSCC, GBM): pERK, ERKSUB, ETV composite. Descriptive."""
import numpy as np, pandas as pd
from lib import *
rows = []
for coh in ['LSCC', 'HNSCC', 'GBM']:
    A = load(coh)
    A['grp_amp'] = np.where(A.EGFR_mut == 1, 'EGFRmut', np.where(A.EGFR_gistic == 2, 'EGFRamp', np.where(A.group == 'RAS', 'RAS', 'WT')))
    print(coh, A.grp_amp.value_counts().to_dict())
    for o in ['stage_RTK', 'stage_ERK', 'pERK_broad', 'stage_ERKSUB', 'ETVc', 'rna_ETV4', 'rna_ETV5', 'rna_DUSP6', 'prot_DUSP6', 'rna_DUSP4']:
        if A[o].notna().sum() < 20: continue
        rows += row(ols(A, o, ['purity', 'prolif'], groups=('EGFRmut', 'EGFRamp'), gcol='grp_amp'), cohort=coh, outcome=o)
T = pd.DataFrame(rows); T.to_csv(OUT + 'tables/T9_other_cohorts.tsv', sep='\t', index=False)
pd.set_option('display.width', 200); print(T[['cohort', 'outcome', 'contrast', 'beta', 'lo', 'hi', 'p', 'n_g', 'n_ref']].round(2).to_string())
