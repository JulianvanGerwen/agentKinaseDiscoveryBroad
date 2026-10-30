"""Shared vectorised 4PL grid fitter (copied from 01_fit.py)."""
import numpy as np
HILLS = np.array([0.3,0.4,0.5,0.6,0.7,0.85,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,5.0])


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


