"""Recompute the headline up/down replication asymmetry with an independent implementation.
Variants: centred vs raw; I19 group set (overlapping, incl. axes) vs disjoint classes; site-weighted vs study-pair-weighted.
Metrics: rep (P same-dir >0.3 in B), base (per-direction base rate in B over co-measured sites), snull (other-class exps of study B, site-matched), lor."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import *
OUT = '/workspace/results/audit_I19/'
X, E = load_exp('inh')
keep = X.notna().sum(1).values >= 6
X = X[keep]
A0 = mask_binders(X, E)
groupsI19 = {c: (E.cls == c).values for c in ['MEK', 'EGFR', 'MET_ALK', 'mTORcat', 'CDK', 'RAF', 'PI3K', 'AKT', 'SRC_ABL', 'DDR']}
groupsI19['AXIS_MAPK'] = (E.cls.isin(['MEK', 'ERK']) | ((E.cls == 'RAF') & E.cell.isin(V600))).values
groupsI19['AXIS_PI3K'] = E.cls.isin(['PI3K', 'PI3K_mTOR', 'mTORcat', 'AKT']).values
groupsI19['AXIS_RTK'] = E.cls.isin(['EGFR', 'MET_ALK', 'RTKmulti', 'FGFR']).values
disjoint = {c: (E.cls == c).values for c in E.cls.unique() if E[E.cls == c].study.nunique() >= 2}
res = {}
for cen in ['centred', 'raw']:
    A = center(A0, E) if cen == 'centred' else A0
    np.save(OUT + f'cache/A_{cen}.npy', A)
    for gname, groups in [('I19groups', groupsI19), ('disjoint', disjoint)]:
        R = pair_rows(A, E, groups)
        R.to_csv(OUT + f'cache/pairs_{cen}_{gname}.tsv.gz', sep='\t', index=False)
        for w in ['sites', 'study_pair']:
            S = summarise(R, weight=w); S.insert(0, 'variant', f'{cen}|{gname}|{w}'); res[(cen, gname, w)] = S
            print(S.round(3).to_string(), flush=True)
pd.concat(res.values()).to_csv(OUT + 'T_headline_recompute.tsv', sep='\t', index=False)
pd.Series(X.index).to_csv(OUT + 'cache/A_sites.tsv', sep='\t', index=False)
