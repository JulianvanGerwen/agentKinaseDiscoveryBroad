"""Task 3 (run with arg "down" for the mirror analysis of decreasing sites): site-level feedback invariance.
Groups: inhibitor classes with >=3 studies + pooled axes (MAPK = MEK/ERK + RAF in BRAF-V600 lines; RTK = EGFR/MET_ALK/RTKmulti/FGFR;
PI3K = PI3K/PI3K_mTOR/mTORcat/AKT).
Values: centred log2FC (within study x cell line where possible). Sites on proteins bound by the drug (primary/Klaeger/Hijazi) are masked.
Per study: mean over the group's experiments; study-level z = mean / (s_site / sqrt(n)), s_site = robust SD of site over all inhibitor experiments (floor 0.25).
'up in study' = mean > 0.5 and z > 2.  Invariant-up: measured in >=3 studies, up in >=3, and up in >=60% of measured studies, none 'down' (mean<-0.5,z<-2).
Context-specific-up: measured in >=3 studies, up in exactly 1 with z>3, and mean<=0.2 in all others.
Null: matched random-condition draws (same per-study counts from other-group experiments of the same study, else global pool), 100 draws -> expected # invariant (FDR)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import *
OUT = "/workspace/results/b_I19/"
rng = np.random.default_rng(7)
E = experiments()
X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]
X = X[X.notna().sum(1) >= 6]
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'gene', 'uniprot'], low_memory=False).set_index('site')
Xc = center(X, E)
gene = S.gene.reindex(Xc.index).fillna('')
A = Xc.values.astype(np.float32)
SIGN = -1 if (len(sys.argv) > 1 and sys.argv[1] == 'down') else 1
SUF = '_down' if SIGN < 0 else ''
A = A * SIGN
# mask own-binder proteins
for j, b in enumerate(E.bset):
    if b: A[gene.isin(b).values, j] = np.nan
med = np.nanmedian(A, 1, keepdims=True)
s = np.maximum(1.4826 * np.nanmedian(np.abs(A - med), 1), 0.25)
v600 = {'COLO205', 'WIDR', 'A375', 'WM239A', 'HT29'}
groups = {c: (E.cls == c).values for c in E.cls.unique() if E[E.cls == c].study.nunique() >= 3}
groups['AXIS_MAPK'] = (E.cls.isin(['MEK', 'ERK']) | ((E.cls == 'RAF') & E.cell_line.isin(v600))).values
groups['AXIS_RTK'] = E.cls.isin(['EGFR', 'MET_ALK', 'RTKmulti', 'FGFR']).values
groups['AXIS_PI3K'] = E.cls.isin(['PI3K', 'PI3K_mTOR', 'mTORcat', 'AKT']).values
studies = E.study.values
def stats_for(idx):
    us = np.unique(studies[idx]); M = np.full((A.shape[0], len(us)), np.nan); Z = M.copy()
    for j, u in enumerate(us):
        sub = A[:, idx & (studies == u)]; n = np.sum(~np.isnan(sub), 1)
        with np.errstate(all='ignore'): m = np.where(n > 0, np.nanmean(sub, 1), np.nan)
        M[:, j] = m; Z[:, j] = m / (s / np.sqrt(np.maximum(n, 1)))
    meas = np.sum(~np.isnan(M), 1); up = (M > 0.5) & (Z > 2); dn = (M < -0.5) & (Z < -2)
    nup = up.sum(1); ndn = dn.sum(1)
    inv = (meas >= 3) & (nup >= 3) & (nup >= 0.6 * meas) & (ndn == 0)
    ctx = (meas >= 3) & (np.sum((M > 0.5) & (Z > 3), 1) == 1) & (np.sum(np.nan_to_num(M, nan=-9) > 0.2, 1) == 1)
    return M, Z, us, meas, nup, ndn, inv, ctx
allrows = []; summ = []
for gname, idx in groups.items():
    M, Z, us, meas, nup, ndn, inv, ctx = stats_for(idx)
    # null
    other = ~idx; ninv_null = []
    for p in range(100):
        pick = np.zeros(len(E), bool)
        for u in np.unique(studies[idx]):
            nK = np.sum(idx & (studies == u)); pool = np.where(other & (studies == u))[0]
            if len(pool) < nK: pool = np.where(other & ~np.isin(studies, np.unique(studies[idx])))[0]
            pick[rng.choice(pool, nK, replace=False)] = True
        ninv_null.append(stats_for(pick)[6].sum())
    nn = np.array(ninv_null)
    summ.append(dict(group=gname, n_exp=int(idx.sum()), n_studies=len(us), n_sites_meas3=int((meas >= 3).sum()), n_invariant_up=int(inv.sum()),
                     null_mean=nn.mean(), null_95=np.percentile(nn, 95), FDR_est=nn.mean() / max(inv.sum(), 1), n_context_up=int(ctx.sum())))
    print(summ[-1], flush=True)
    with np.errstate(all='ignore'): mmean = np.nanmean(M, 1)
    for i in np.where(inv | ctx)[0]:
        allrows.append(dict(group=gname, site=Xc.index[i], gene=gene.iloc[i], kind='invariant' if inv[i] else 'context', meas_studies=int(meas[i]), up_studies=int(nup[i]),
                            mean_study_lfc=mmean[i], per_study=';'.join(f'{u}:{m:.2f}' for u, m in zip(us, M[i]) if not np.isnan(m))))
    # store per-site stats for predictors (sites measured in >=3 studies)
    keep = meas >= 3
    pd.DataFrame({'site': Xc.index[keep], 'meas': meas[keep], 'nup': nup[keep], 'ndn': ndn[keep], 'inv': inv[keep], 'ctx': ctx[keep], 'mean': mmean[keep]}).to_parquet(CA + f'site_stats{SUF}_{gname}.parquet')
pd.DataFrame(summ).to_csv(OUT + f'tables/T3_invariance_summary{SUF}.tsv', sep='\t', index=False)
pd.DataFrame(allrows).to_csv(OUT + f'tables/T3_invariant_and_context_sites{SUF}.tsv', sep='\t', index=False)
