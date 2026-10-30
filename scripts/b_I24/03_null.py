"""Random-site null: replace output stage with residue/size-matched random site sets (same cohort);
also replace input stage. Compare observed partial r (adj purity, prolif) with null distribution."""
from common import *
rng = np.random.default_rng(1)
import sys
NPC = int(sys.argv[1]) if len(sys.argv) > 1 else 0
A = load_stages('tc_padj'); NR = 200
PCS = pd.read_parquet(f'{OUT}/tables/phospho_pcs.parquet')
stage_prot = set(U.values())
def pr(x, y, C):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30: return np.nan
    Cc = np.c_[np.ones(ok.sum()), C[ok]]
    rx = x[ok] - Cc @ np.linalg.lstsq(Cc, x[ok], rcond=None)[0]; ry = y[ok] - Cc @ np.linalg.lstsq(Cc, y[ok], rcond=None)[0]
    return np.corrcoef(rx, ry)[0, 1]
rows = []
for c in COHORTS:
    D = A[A.cohort == c]; X = cl.matrix(c, 'input_tc_padj')[D.index]
    X = X[[k.split('_')[0] not in stage_prot for k in X.index]]
    cov = X.notna().mean(1); X = X[cov >= 0.5]
    res = np.array([k.split('_')[1][0] for k in X.index]); C = np.c_[D[['purity', 'prolif']].values, PCS.loc[D.index].values[:, :NPC]]
    Xz = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
    for e in EDGES:
        xs, ys = e
        if D[xs].notna().sum() < 30 or D[ys].notna().sum() < 30: continue
        obs = pr(D[xs].values, D[ys].values, C)
        if ys == 'ERKSUB': nsite, comp = 25, None
        else:
            sites = [s for s in STAGES[ys] if s in cl.matrix(c, 'input_tc_padj').index]
            comp = [s.split('_')[1][0] for s in sites]; nsite = len(sites)
        null = []
        for _ in range(NR):
            if comp is None:
                pick = rng.choice(len(X), nsite, replace=False); v = Xz.iloc[pick].mean(0).values
            else:
                pick = [rng.choice(np.where(res == r)[0]) for r in comp]; v = X.iloc[pick].mean(0).values
            null.append(pr(D[xs].values, v, C))
        null = np.array(null)
        rows.append(dict(cohort=c, edge=ename(e), obs_r=obs, null_mean=np.nanmean(null), null_sd=np.nanstd(null),
                         z=(obs - np.nanmean(null)) / np.nanstd(null), pct=np.nanmean(null < obs)))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/tables/T3_random_site_null_pc{NPC}.tsv', sep='\t', index=False)
S_ = R.groupby('edge').agg(obs_r=('obs_r', 'median'), null_mean=('null_mean', 'median'), z_median=('z', 'median'),
                           n_cohort_pct95=('pct', lambda p: int((p > 0.95).sum())), k=('pct', 'size'))
print(S_.round(3).to_string())
