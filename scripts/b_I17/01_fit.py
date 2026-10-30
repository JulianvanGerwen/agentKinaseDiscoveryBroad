"""Vectorised 4-parameter log-logistic fits per site per dose experiment, with permuted-dose null.
Model: y = front + delta * s(x),  s = 1/(1+10^(h*(logEC50 - log10 x)))  (delta<0 = down-regulated).
For each (logEC50, h) on a grid the model is linear in (front, delta) -> closed-form weighted LS
over all sites at once (NaN-masked); best grid point per site = 4PL fit (grid step 0.05 log, 16 Hill values).
Statistic: F = ((SS0-SSE)/3)/(SSE/(n-4)) vs flat model. Null: same fit after jointly permuting
dose labels (columns) within experiment (preserves site-site correlation), NPERM permutations.
Output: cache/b_I17/fits.parquet, cache/b_I17/null_F.parquet
"""
import os, sys, time
for v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]: os.environ[v] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load

OUT = "/workspace/cache/b_I17"
NPERM = int(os.environ.get("NPERM", 20))
HILLS = np.array([0.3,0.4,0.5,0.6,0.7,0.85,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,5.0])
E = pd.read_csv(f"{OUT}/experiments.tsv", sep="\t")
E = E[(E.pclass != "time_course") & (E.n_doses >= 6)].reset_index(drop=True)
only = os.environ.get("ONLY")
if only: E = E[E.exp_id.isin(only.split(","))]


def fit_matrix(Y, lx, lo, hi):
    """Y: sites x doses (NaN allowed); lx: log10 dose per column. Returns dict of arrays."""
    M = np.isfinite(Y).astype(np.float64); Y0 = np.where(M > 0, Y, 0.0)
    n = M.sum(1); sy = Y0.sum(1); syy = (Y0 ** 2).sum(1)
    ss0 = syy - sy ** 2 / np.maximum(n, 1)
    grid_ec = np.arange(lo - 0.5, hi + 0.5 + 1e-9, 0.05)
    EC, H = np.meshgrid(grid_ec, HILLS, indexing="ij"); EC = EC.ravel(); H = H.ravel()
    S = 1.0 / (1.0 + 10 ** (H[None, :] * (EC[None, :] - lx[:, None])))  # doses x G
    ss_ = M @ S; ss2 = M @ (S ** 2); ssy = Y0 @ S                     # sites x G
    det = n[:, None] * ss2 - ss_ ** 2
    det = np.where(np.abs(det) < 1e-9, np.nan, det)
    d = (n[:, None] * ssy - ss_ * sy[:, None]) / det
    f = (sy[:, None] - d * ss_) / n[:, None]
    sse = syy[:, None] - f * sy[:, None] - d * ssy
    sse = np.where(np.isfinite(sse), sse, np.inf)
    j = np.argmin(sse, 1); r = np.arange(len(j))
    inr = (EC >= lo) & (EC <= hi)
    sse_h1 = np.min(np.where((H == 1.0)[None, :] & inr[None, :], sse, np.inf), 1)  # Hill fixed at 1
    sse_inr = np.min(np.where(inr[None, :], sse, np.inf), 1)
    best = sse[r, j]
    F = ((ss0 - best) / 3) / (np.maximum(best, 1e-12) / np.maximum(n - 4, 1))
    return dict(n=n, logEC50=EC[j], hill=H[j], front=f[r, j], delta=d[r, j], sse=best, ss0=ss0,
                R2=1 - best / np.maximum(ss0, 1e-12), F=F, sse_h1=sse_h1, sse_inr=sse_inr)


rng = np.random.default_rng(17)
allfits, nulls = [], []
X = None
for _, e in E.iterrows():
    t0 = time.time()
    conds = e.conditions.split(";"); doses = np.array([float(v) for v in e.doses.split(";")])
    Xe = load.sites_matrix(conds)
    Xe = Xe[Xe.notna().sum(axis=1) >= 6]
    nz = doses[doses > 0]; lo, hi = np.log10(nz.min()), np.log10(nz.max())
    lx = np.where(doses > 0, np.log10(np.maximum(doses, 1e-12)), lo - 1.0)
    Y = Xe.values.astype(np.float64)
    res = fit_matrix(Y, lx, lo, hi)
    df = pd.DataFrame(res); df.insert(0, "site", Xe.index.values); df.insert(0, "exp_id", e.exp_id)
    df["ec_in_range"] = (df.logEC50 >= lo) & (df.logEC50 <= hi)
    df["max_abs"] = np.nanmax(np.abs(Y), 1)
    allfits.append(df)
    for p in range(NPERM):
        perm = rng.permutation(len(lx))
        rp = fit_matrix(Y, lx[perm], lo, hi)
        inr = (rp["logEC50"] >= lo) & (rp["logEC50"] <= hi)
        nulls.append(pd.DataFrame(dict(exp_id=e.exp_id, perm=p, F=rp["F"].astype(np.float32), delta=rp["delta"].astype(np.float32),
                                       R2=rp["R2"].astype(np.float32), ec_in_range=inr)))
    print(f"{e.exp_id}: {len(Xe)} sites, {len(conds)} cond, {time.time()-t0:.1f}s", flush=True)

pd.concat(allfits).to_parquet(f"{OUT}/fits{'_'+only if only else ''}.parquet")
pd.concat(nulls).to_parquet(f"{OUT}/null_F{'_'+only if only else ''}.parquet")
