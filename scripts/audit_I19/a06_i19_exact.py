"""Re-implementation of the I19 study-level metric from the audit build (not from b_I19 cache):
per group & study: mean over the group's experiments (centred, binders masked); discovery |mean_A| in bin & |z|>2; replication same-sign mean_B > 0.3.
Null: mean over an equal-size random draw of other-group experiments in study B (I19: 1 draw; here 20 draws averaged).
Reports (i) I19-style aggregation (rep over all pairs, null over the subset with a null) and (ii) matched aggregation (only pairs with a null)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import *
OUT = '/workspace/results/audit_I19/'
X, E = load_exp('inh'); X = X[X.notna().sum(axis=1).values >= 6]
A = np.load(OUT + 'cache/A_centred.npy')
s = np.maximum(1.4826 * np.nanmedian(np.abs(A - np.nanmedian(A, 1, keepdims=True)), 1), 0.25)
groups = {c: (E.cls == c).values for c in ['MEK', 'EGFR', 'MET_ALK', 'mTORcat', 'CDK', 'RAF', 'PI3K', 'AKT', 'SRC_ABL', 'DDR']}
groups['AXIS_MAPK'] = (E.cls.isin(['MEK', 'ERK']) | ((E.cls == 'RAF') & E.cell.isin(V600))).values
groups['AXIS_PI3K'] = E.cls.isin(['PI3K', 'PI3K_mTOR', 'mTORcat', 'AKT']).values
groups['AXIS_RTK'] = E.cls.isin(['EGFR', 'MET_ALK', 'RTKmulti', 'FGFR']).values
studies = E.study.values; rng = np.random.default_rng(11); rows = []
for g, idx in groups.items():
    us = np.unique(studies[idx]); M, Z, N = {}, {}, {}
    for u in us:
        m_ = idx & (studies == u); sub = A[:, m_]; n = (~np.isnan(sub)).sum(1)
        with np.errstate(all='ignore'): M[u] = np.where(n > 0, np.nanmean(sub, 1), np.nan)
        Z[u] = M[u] / (s / np.sqrt(np.maximum(n, 1)))
        pool = np.where(~idx & (studies == u))[0]
        if len(pool) >= m_.sum():
            draws = []
            for _ in range(20):
                pk = rng.choice(pool, m_.sum(), replace=False)
                with np.errstate(all='ignore'): draws.append(np.nanmean(A[:, pk], 1))
            N[u] = np.stack(draws, 1)
    for a in us:
        for b in us:
            if a == b: continue
            both = ~np.isnan(M[a]) & ~np.isnan(M[b])
            for sg, lab in [(1, 'up'), (-1, 'down')]:
                ma, mb, za = sg * M[a], sg * M[b], sg * Z[a]
                for lo, hi in BINS:
                    d = both & (ma > lo) & (ma <= hi) & (za > 2)
                    if d.sum() < 5: continue
                    r = dict(group=g, si=a, sj=b, dir=lab, bin=f'{lo}-{hi}', n_disc=int(d.sum()), rep=np.mean(mb[d] > 0.3), base=np.mean(mb[both] > 0.3))
                    if b in N:
                        nb = sg * N[b][d]; r['null'] = np.nanmean(np.where(np.isnan(nb), np.nan, nb > 0.3))
                    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(OUT + 'cache/pairs_i19exact.tsv.gz', sep='\t', index=False)
wm = lambda d, c: np.average(d[c].dropna(), weights=d.loc[d[c].notna(), 'n_disc'])
out = []
for (dr, b), d in R.groupby(['dir', 'bin']):
    dn = d[d.null.notna()]; dd = d[~d.group.str.startswith('AXIS')]; ddn = dd[dd.null.notna()]
    out.append(dict(dir=dr, bin=b, rep_I19style=wm(d, 'rep'), null_I19style=wm(d, 'null'), rep_matched=wm(dn, 'rep'), null_matched=wm(dn, 'null'), base_matched=wm(dn, 'base'),
                    rep_disjoint_matched=wm(ddn, 'rep'), null_disjoint_matched=wm(ddn, 'null'), frac_pairs_with_null=d.null.notna().mean()))
out = pd.DataFrame(out); print(out.round(3).to_string()); out.to_csv(OUT + 'T_i19exact_recompute.tsv', sep='\t', index=False)
