import numpy as np, pandas as pd
from scipy import stats

CACHE = '/workspace/results/b_I20/cache'
RES = '/workspace/results/b_I20'
COHORTS = ['BRCA', 'CCRCC', 'COAD', 'GBM', 'HNSCC', 'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC']


def load(C, name):
    return pd.read_parquet(f'{CACHE}/{C}/{name}.parquet')


def resid_cov(Z, cov):
    """Z: (n x k) complete; residualise on [1, cov]"""
    X = np.column_stack([np.ones(len(cov)), cov])
    B = np.linalg.lstsq(X, Z, rcond=None)[0]
    return Z - X @ B


def rowwise_ols_t(Y, z, cov):
    """per row g: Y[g] ~ 1 + cov + z[g] (masked by NaN in Y[g] or z[g]). Y,z: (G x n). returns beta,t,n"""
    G, n = Y.shape
    M = np.isfinite(Y) & np.isfinite(z)
    X = np.concatenate([np.ones((G, n, 1)), np.broadcast_to(cov[None], (G, n, cov.shape[1])), np.nan_to_num(z)[:, :, None]], 2)
    p = X.shape[2]
    Xm = X * M[:, :, None]
    Y0 = np.where(M, Y, 0)
    XtX = np.einsum('gnj,gnk->gjk', Xm, X) + 1e-9 * np.eye(p)
    Xty = np.einsum('gnj,gn->gj', Xm, Y0)
    inv = np.linalg.inv(XtX)
    B = np.einsum('gjk,gk->gj', inv, Xty)
    R = np.where(M, Y0 - np.einsum('gnj,gj->gn', X, B), 0)
    nn = M.sum(1)
    s2 = (R ** 2).sum(1) / np.maximum(nn - p, 1)
    se = np.sqrt(s2 * inv[:, -1, -1])
    return B[:, -1], B[:, -1] / se, nn


def masked_t(Z, Y):
    """Z: (n x K) complete regressors (already cov-residualised); Y: (n x S) with NaN.
    simple regression of y_s on z_k using samples where y_s observed (centred within mask).
    returns beta (K x S), t (K x S), n (S)"""
    M = np.isfinite(Y).astype(np.float64)
    Y0 = np.nan_to_num(Y)
    n = M.sum(0)
    Sz = Z.T @ M; Szz = (Z ** 2).T @ M; Szy = Z.T @ Y0
    Sy = Y0.sum(0); Syy = (Y0 ** 2).sum(0)
    cxx = Szz - Sz ** 2 / n
    cxy = Szy - Sz * Sy / n
    cyy = Syy - Sy ** 2 / n
    beta = cxy / cxx
    r = cxy / np.sqrt(cxx * cyy)
    r = np.clip(r, -0.999999, 0.999999)
    t = r * np.sqrt((n - 3) / (1 - r ** 2))   # df n-2 minus 1 for cov already removed (approx)
    return beta, t, n


def auroc(scores, labels):
    labels = np.asarray(labels, bool); scores = np.asarray(scores, float)
    ok = np.isfinite(scores); scores, labels = scores[ok], labels[ok]
    npos = labels.sum(); nneg = (~labels).sum()
    if npos == 0 or nneg == 0:
        return np.nan
    rk = stats.rankdata(scores)
    return (rk[labels].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = np.isfinite(p)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        q[ok] = multipletests(p[ok], method='fdr_bh')[1]
    return q
