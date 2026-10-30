"""Task 3b/4: up-vs-down asymmetry in cross-study replication, magnitude-matched.
For each group and each ordered pair of studies (A discovery, B replication), for sites measured in both:
P(site responds in B in the same direction | it responds in A with |study mean| in bin), for up and down separately.
Response threshold in B: same-sign study mean > 0.3 (lenient).  Discovery bins on |mean_A|: 0.5-1, 1-2, >2 (and z>2).
Null: same computation with discovery in A and the replication values drawn from a matched random set of other-group experiments in B."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import *
OUT = "/workspace/results/b_I19/"
E = experiments()
LOWDOWN = len(sys.argv) > 1 and sys.argv[1] == 'lowdown'
if LOWDOWN:  # artefact control: drop experiments in the top quartile of fraction strongly-down sites
    G = pd.read_csv(CA + 'exp_globalstats.tsv', sep='\t', index_col=0)
    E = E[E.exp.map(G.frac_down1) <= G.frac_down1.quantile(0.75)].reset_index(drop=True)
SUF = '_lowdown' if LOWDOWN else ''
RAW = len(sys.argv) > 1 and sys.argv[1] == 'raw'
if RAW:  # control for centring artefacts: uncentred vs-DMSO values; drop study 32034161 (shared EGF offset)
    E = E[E.study != '32034161'].reset_index(drop=True); E['centerable'] = False; SUF = '_raw'
X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]; X = X[X.notna().sum(1) >= 6]
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'gene'], low_memory=False).set_index('site')
Xc = center(X, E); gene = S.gene.reindex(Xc.index).fillna('')
A = Xc.values.astype(np.float32)
for j, b in enumerate(E.bset):
    if b: A[gene.isin(b).values, j] = np.nan
s = np.maximum(1.4826 * np.nanmedian(np.abs(A - np.nanmedian(A, 1, keepdims=True)), 1), 0.25)
v600 = {'COLO205', 'WIDR', 'A375', 'WM239A', 'HT29'}
groups = {c: (E.cls == c).values for c in ['MEK', 'EGFR', 'MET_ALK', 'mTORcat', 'CDK', 'RAF', 'PI3K', 'AKT', 'SRC_ABL', 'DDR']}
groups['AXIS_MAPK'] = (E.cls.isin(['MEK', 'ERK']) | ((E.cls == 'RAF') & E.cell_line.isin(v600))).values
groups['AXIS_PI3K'] = E.cls.isin(['PI3K', 'PI3K_mTOR', 'mTORcat', 'AKT']).values
groups['AXIS_RTK'] = E.cls.isin(['EGFR', 'MET_ALK', 'RTKmulti', 'FGFR']).values
studies = E.study.values
bins = [(0.5, 1), (1, 2), (2, 99)]
rows = []
rng = np.random.default_rng(3)
for g, idx in groups.items():
    us = np.unique(studies[idx]); M = {}; Zs = {}; Nul = {}
    for u in us:
        m_ = idx & (studies == u); sub = A[:, m_]; n = np.sum(~np.isnan(sub), 1)
        with np.errstate(all='ignore'): M[u] = np.where(n > 0, np.nanmean(sub, 1), np.nan)
        Zs[u] = M[u] / (s / np.sqrt(np.maximum(n, 1)))
        pool = np.where(~idx & (studies == u))[0]
        if len(pool) >= m_.sum():
            pk = rng.choice(pool, m_.sum(), replace=False)
            with np.errstate(all='ignore'): Nul[u] = np.nanmean(A[:, pk], 1)
    for a in us:
        for b in us:
            if a == b: continue
            both = ~np.isnan(M[a]) & ~np.isnan(M[b])
            for sign, lab in [(1, 'up'), (-1, 'down')]:
                ma = sign * M[a]; mb = sign * M[b]; za = sign * Zs[a]
                for lo, hi in bins:
                    d = both & (ma > lo) & (ma <= hi) & (za > 2)
                    if d.sum() < 5: continue
                    r = dict(group=g, A=a, B=b, dir=lab, bin=f'{lo}-{hi}', n_disc=int(d.sum()), rep=float(np.mean(mb[d] > 0.3)), mean_B=float(np.mean(mb[d])))
                    if b in Nul:
                        nb = sign * Nul[b]; dd = d & ~np.isnan(nb); r['rep_null'] = float(np.mean(nb[dd] > 0.3)) if dd.sum() else np.nan
                    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(OUT + f'tables/T3_replication_pairs{SUF}.tsv', sep='\t', index=False)
# summarise: weighted mean replication over study pairs
def wm(d, c): return np.average(d[c].dropna(), weights=d.loc[d[c].notna(), 'n_disc']) if d[c].notna().any() else np.nan
Sm = R.groupby(['group', 'dir', 'bin']).apply(lambda d: pd.Series({'pairs': len(d), 'n_disc': d.n_disc.sum(), 'rep': wm(d, 'rep'), 'rep_null': wm(d, 'rep_null'), 'mean_B': wm(d, 'mean_B')})).reset_index()
Sm['rep_excess'] = Sm.rep - Sm.rep_null
Sm.to_csv(OUT + f'tables/T3_replication_summary{SUF}.tsv', sep='\t', index=False)
pd.set_option('display.width', 200); pd.set_option('display.max_rows', 200)
print(Sm.pivot_table(index=['group', 'bin'], columns='dir', values=['rep', 'rep_null', 'n_disc']).round(3).to_string())
# pooled over groups
P = R.groupby(['dir', 'bin']).apply(lambda d: pd.Series({'n_disc': d.n_disc.sum(), 'rep': wm(d, 'rep'), 'rep_null': wm(d, 'rep_null')})).reset_index(); print(P)
