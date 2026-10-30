"""Task 3c: which context axis breaks 'up' invariance?  Experiment-pair replication by pair type:
  same_drug_diff_cell (same study), same_class_diff_drug_same_cell (same study), same_class_diff_study,
  and control pairs: diff_class_same_cell (same study) = null.
Discovery: centred value > 1 (up) or < -1 (down) in experiment 1; replication: same-sign value > 0.5 in experiment 2.
Own-binder proteins masked."""
import sys, numpy as np, pandas as pd, itertools
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import *
OUT = "/workspace/results/b_I19/"
E = experiments().reset_index(drop=True)
X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]; X = X[X.notna().sum(1) >= 6]
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'gene'], low_memory=False).set_index('site')
Xc = center(X, E); gene = S.gene.reindex(Xc.index).fillna('')
A = Xc.values.astype(np.float32)
for j, b in enumerate(E.bset):
    if b: A[gene.isin(b).values, j] = np.nan
skip = {'other', 'RTKmulti'}
rows = []
for i, j in itertools.permutations(range(len(E)), 2):
    ei, ej = E.loc[i], E.loc[j]
    if ei.cls in skip or ej.cls in skip: continue
    if ei.study == ej.study:
        if ei.agent == ej.agent and ei.cell_line != ej.cell_line: t = 'same_drug_diff_cell'
        elif ei.cls == ej.cls and ei.agent != ej.agent and ei.cell_line == ej.cell_line: t = 'same_class_diff_drug_same_cell'
        elif ei.cls == ej.cls and ei.agent != ej.agent and ei.cell_line != ej.cell_line: t = 'same_class_diff_drug_diff_cell'
        elif ei.cls != ej.cls and ei.cell_line == ej.cell_line: t = 'NULL_diff_class_same_cell'
        elif ei.cls != ej.cls: t = 'NULL_diff_class_diff_cell'
        else: continue
    else:
        if ei.cls == ej.cls: t = 'same_class_diff_study'
        else: continue  # cross-study diff-class null is huge; skip
    a, b = A[:, i], A[:, j]; m = ~np.isnan(a) & ~np.isnan(b)
    if m.sum() < 500: continue
    a, b = a[m], b[m]
    up = a > 1; dn = a < -1
    rows.append(dict(type=t, e1=ei.exp, e2=ej.exp, cls=ei.cls, n_up=int(up.sum()), n_dn=int(dn.sum()),
                     rep_up=float(np.mean(b[up] > 0.5)) if up.sum() >= 5 else np.nan, rep_dn=float(np.mean(b[dn] < -0.5)) if dn.sum() >= 5 else np.nan,
                     bg_up=float(np.mean(b > 0.5)), bg_dn=float(np.mean(b < -0.5))))
R = pd.DataFrame(rows); R.to_csv(OUT + 'tables/T3_pairtype_replication.tsv', sep='\t', index=False)
R['ex_up'] = R.rep_up - R.bg_up; R['ex_dn'] = R.rep_dn - R.bg_dn
g = R.groupby('type').agg(pairs=('e1', 'size'), rep_up=('rep_up', 'median'), rep_dn=('rep_dn', 'median'), ex_up=('ex_up', 'median'), ex_dn=('ex_dn', 'median'),
                          bg_up=('bg_up', 'median'), bg_dn=('bg_dn', 'median'))
print(g.round(3).to_string())
g.to_csv(OUT + 'tables/T3_pairtype_summary.tsv', sep='\t')
# per-class for same_drug_diff_cell and same_class_diff_drug_same_cell
print(R[R.type.isin(['same_drug_diff_cell', 'same_class_diff_drug_same_cell', 'same_class_diff_study'])].groupby(['type', 'cls'])[['ex_up', 'ex_dn']].median().round(3).to_string())
