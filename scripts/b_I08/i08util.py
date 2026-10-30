import numpy as np, pandas as pd
from scipy import stats
import sys
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
R = cl.ROOT; OUT = '/workspace/results/b_I08'


def rank(x):
    x = np.asarray(x, float); r = np.full_like(x, np.nan); m = np.isfinite(x)
    r[m] = stats.rankdata(x[m]); return r


def resid(y, C):
    """residualise y on C (n x k, complete) using finite y rows; returns array with NaN where y NaN."""
    m = np.isfinite(y); out = np.full_like(y, np.nan)
    if m.sum() < C.shape[1] + 5: return out
    X = np.column_stack([np.ones(m.sum()), C[m]])
    b, *_ = np.linalg.lstsq(X, y[m], rcond=None); out[m] = y[m] - X @ b; return out


def pcor(a, b):
    m = np.isfinite(a) & np.isfinite(b); n = int(m.sum())
    if n < 15: return np.nan, n
    x, y = a[m], b[m]
    if x.std() == 0 or y.std() == 0: return np.nan, n
    return float(np.corrcoef(x, y)[0, 1]), n


def meta(r, n, k=0):
    """DerSimonian-Laird random-effects meta of correlations via Fisher z. k = #covariates partialled.
    Returns pooled r, lo, hi, p, I2, n_cohorts"""
    r = np.asarray(r, float); n = np.asarray(n, float); m = np.isfinite(r) & (n - 3 - k > 0)
    r, n = np.clip(r[m], -0.999, 0.999), n[m]
    if len(r) < 2: return dict(r=np.nan, lo=np.nan, hi=np.nan, p=np.nan, I2=np.nan, k=len(r))
    z = np.arctanh(r); v = 1 / (n - 3 - k); w = 1 / v
    zf = (w * z).sum() / w.sum(); Q = (w * (z - zf) ** 2).sum(); df = len(z) - 1
    tau2 = max(0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum()))
    ws = 1 / (v + tau2); zr = (ws * z).sum() / ws.sum(); se = np.sqrt(1 / ws.sum())
    p = 2 * stats.norm.sf(abs(zr / se)); I2 = max(0, (Q - df) / Q) if Q > 0 else 0
    return dict(r=np.tanh(zr), lo=np.tanh(zr - 1.96 * se), hi=np.tanh(zr + 1.96 * se), p=p, I2=I2, k=len(r))


def bh(p):
    p = np.asarray(p, float); out = np.full_like(p, np.nan); m = np.isfinite(p)
    if m.sum(): out[m] = stats.false_discovery_control(p[m])
    return out


def covariates(c, T):
    """purity (WES, median-imputed), ESTIMATE score, proliferation (mean z of E2F + G2M ssGSEA) for tumour ids T."""
    cases = [t[:-2] for t in T]
    clin = cl.table(c, 'clinical').drop_duplicates('case_id').set_index('case_id').reindex(cases)
    ann = cl.table(c, 'annotations').drop_duplicates('case_id').set_index('case_id').reindex(cases)
    pur = clin.WES_purity.astype(float)
    est = ann.ESTIMATE_ESTIMATEScore.astype(float)
    pro = ann[['HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT']].astype(float)
    pro = ((pro - pro.mean()) / pro.std()).mean(1)
    C = pd.DataFrame({'purity': pur.fillna(pur.median()).values, 'estimate': est.fillna(est.median()).values,
                      'prolif': pro.fillna(pro.median()).values}, index=T)
    return C.fillna(0)
