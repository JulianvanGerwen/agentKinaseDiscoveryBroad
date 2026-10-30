"""Leave-one-class-out and leave-one-study-out for dLOR(down-up) in inhibitors (primary) and activating stimulations."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import *
OUT = '/workspace/results/audit_I19/'
Ri = pd.read_csv(OUT + 'cache/pairs_raw_disjoint.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
Rs = pd.read_csv(OUT + 'cache/pairs_stim.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
Rs = Rs[Rs.group.isin(['GF_RTK', 'PPase_inh', 'DNA_damage', 'GPCR', 'TNF', 'PMA'])]
def dl(R, b):
    sp = R[R.bin == b].groupby(['si', 'sj', 'dir']).lor.mean().unstack('dir').dropna(); return (sp['down'] - sp['up']).mean()
rows = []
for lab, R in [('inhibitor', Ri), ('stim_activating', Rs)]:
    for b in ['1-2', '2-99']:
        full = dl(R, b)
        lc = {c: dl(R[R.group != c], b) for c in R.group.unique()}
        ls = {s: dl(R[(R.si != s) & (R.sj != s)], b) for s in np.unique(np.r_[R.si, R.sj])}
        rows.append(dict(set=lab, bin=b, full=full, LOCO_min=min(lc.values()), LOCO_min_drop=min(lc, key=lc.get), LOCO_max=max(lc.values()),
                         LOSO_min=min(ls.values()), LOSO_min_drop=min(ls, key=ls.get), LOSO_max=max(ls.values()), LOSO_max_drop=max(ls, key=ls.get)))
        if lab == 'inhibitor':
            for drop in [['MEK'], ['MEK', 'ERK', 'RAF']]:
                rows.append(dict(set=f'inhibitor_minus_{"+".join(drop)}', bin=b, full=dl(R[~R.group.isin(drop)], b)))
out = pd.DataFrame(rows); print(out.round(3).to_string()); out.to_csv(OUT + 'T_loo.tsv', sep='\t', index=False)
# bootstrap CI without MEK
for b in ['1-2', '2-99']:
    print('no MEK', b, boot_diff(Ri[Ri.group != 'MEK'], 'lor', b))
