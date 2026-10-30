"""H4 normalisation axis, site-level: cross-study replication after per-experiment median centring / robust scaling."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I25'); from common_i25 import *
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import pair_rows
A, E = load_test('inh'); Av = A.values.astype(np.float32)
groups = {c: (E.cls == c).values for c in E.cls.unique() if E[E.cls == c].study.nunique() >= 2}
med = np.nanmedian(Av, 0); mad = np.nanmedian(np.abs(Av - med), 0) * 1.4826
rows = []
for lab, M in [('raw', Av), ('median_centred', Av - med), ('robust_z_scaled_to_median_mad', (Av - med) / mad * np.nanmedian(mad))]:
    R = pair_rows(M, E, groups, bins=[(1, 2), (1, 99)], null_pool=None)
    sp = R.groupby(['bin', 'dir', 'si', 'sj']).lor.mean().groupby(['bin', 'dir']).mean().unstack('dir')
    for b in sp.index: rows.append(dict(version=lab, bin=b, up=sp.loc[b, 'up'], down=sp.loc[b, 'down'], gap=sp.loc[b, 'down'] - sp.loc[b, 'up']))
    print(lab, flush=True)
T = pd.DataFrame(rows); T.to_csv('/workspace/results/b_I25/tables/T2b_H4_normalisation.tsv', sep='\t', index=False); print(T.round(3))
