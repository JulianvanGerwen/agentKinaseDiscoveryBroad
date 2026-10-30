"""Audit of I24 candidate E (independent code). Stage builder + IV decomposition utilities."""
import warnings; warnings.filterwarnings("ignore")
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
OUT = '/workspace/results/audit_I24'
COH = cl.COHORTS
ACC = dict(EGFR='P00533', ERBB2='P04626', MET='P08581', SHC1='P29353', GAB1='Q13480', PTPN11='Q06124', MAPK3='P27361', MAPK1='P28482')
def S(g, *s): return [f'{ACC[g]}_{x}' for x in s]
SITES = {
    'RTK': S('EGFR', 'Y1092', 'Y1172', 'Y1197', 'Y1110', 'Y1016') + S('ERBB2', 'Y1248', 'Y877') + S('MET', 'Y1234', 'Y1235'),
    'EGFRpY': S('EGFR', 'Y1092', 'Y1172', 'Y1197', 'Y1110', 'Y1016'),
    'ERBB2pY': S('ERBB2', 'Y1248', 'Y877'),
    'ADAPT': S('SHC1', 'Y349', 'Y350', 'Y427') + S('GAB1', 'Y627', 'Y659') + S('PTPN11', 'Y542', 'Y580'),
    'ADAPT_SHC1': S('SHC1', 'Y349', 'Y350', 'Y427'),
    'ADAPT_PG': S('PTPN11', 'Y542', 'Y580') + S('GAB1', 'Y627', 'Y659'),
    'ERK': S('MAPK3', 'T202', 'Y204') + S('MAPK1', 'T185', 'Y187'),
    'ERK3': S('MAPK3', 'T202', 'Y204'), 'ERK2': S('MAPK1', 'T185', 'Y187'),
}
# SFK readouts (activation loop of SRC/YES/FYN not quantified; use classic SFK substrate sites)
SFK_SUB = ['Q9H5V8_Y734', 'Q9H5V8_Y743', 'Q9H5V8_Y762',  # CDCP1
           'Q05397_Y576', 'Q05397_Y577', 'Q05397_Y925', 'Q05397_Y861',  # PTK2 (FAK) SFK sites
           'P49023_Y118', 'P49023_Y31',  # PXN
           'Q03135_Y14',  # CAV1
           'Q14247_Y421', 'Q14247_Y470', 'Q14247_Y486',  # CTTN
           'P56945_Y249', 'P56945_Y410', 'P56945_Y387']  # BCAR1
PATH_PROT = set(ACC.values())


def zs(x):
    x = np.asarray(x, float); return (x - np.nanmean(x)) / np.nanstd(x)


def composite(X, sites):
    ss = [s for s in sites if s in X.index]
    if not ss: return pd.Series(np.nan, index=X.columns), []
    return X.loc[ss].mean(0, skipna=True), ss


def covariates(c, idx):
    case = [s[:-2] for s in idx]
    clin = cl.clinical().drop_duplicates('case_id').set_index('case_id')
    ann = cl.table(c, 'annotations').drop_duplicates('case_id').set_index('case_id')
    D = pd.DataFrame(index=idx)
    D['purity'] = clin.reindex(case)['WES_purity'].values
    D['prolif'] = ann.reindex(case)['HALLMARK_E2F_TARGETS'].values
    for k in ['ESTIMATE_ImmuneScore', 'ESTIMATE_StromalScore']:
        D[k.split('_')[1].replace('Score', '')] = ann.reindex(case)[k].values if k in ann else np.nan
    for v in D.columns: D[v] = D[v].fillna(D[v].median())
    D['case_id'] = case
    return D


def build(c, variant):
    T = cl.tumour_ids(c); X = cl.matrix(c, f'input_{variant}'); X = X[[t for t in T if t in X.columns]]
    D = covariates(c, X.columns); used = {}
    for k, s in SITES.items():
        D[k], used[k] = composite(X, s)
    # background pY: non-pathway pY sites with >=50% coverage, z-scored, mean
    py = [k for k in X.index if k.split('_')[1][0] == 'Y' and k.split('_')[0] not in PATH_PROT]
    Y = X.loc[py]; Y = Y[Y.notna().mean(1) >= 0.5]; Yz = Y.sub(Y.mean(1), axis=0).div(Y.std(1), axis=0)
    D['bgPY'] = Yz.mean(0); used['bgPY'] = [f'{len(Y)} sites']
    sfk = [s for s in SFK_SUB if s in X.index]
    Z = X.loc[sfk]; Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1), axis=0)
    D['SFK'] = Z.mean(0).where(Z.notna().sum(0) >= 2); used['SFK'] = sfk
    return D, used


def _ols(y, X):
    return np.linalg.lstsq(X, y, rcond=None)[0]


def decomp(D, x='RTK', m='ADAPT', y='ERK', cov=('purity', 'prolif'), nboot=1000, seed=0):
    """Standardized: FS: m~x+cov; RF: y~x+cov; IV=RF/FS (=coef of RTK-driven part, I24 a_rtk);
    OTHER: coef of m in y~m+x+cov (=coef of RTK-independent residual part, I24 a_oth); total y~m+cov; x->m, x->y."""
    cols = [x, m, y] + list(cov)
    d = D[cols].dropna()
    n = len(d)
    if n < 30 or d[x].std() == 0: return None
    Z = d.apply(zs).values
    xi, mi, yi = Z[:, 0], Z[:, 1], Z[:, 2]; C = np.c_[np.ones(n), Z[:, 3:]]
    def est(idx):
        xx, mm, yy, CC = xi[idx], mi[idx], yi[idx], C[idx]
        fs = _ols(mm, np.c_[xx, CC])[0]; rf = _ols(yy, np.c_[xx, CC])[0]
        oth = _ols(yy, np.c_[mm, xx, CC])[0]; tot = _ols(yy, np.c_[mm, CC])[0]
        return np.array([fs, rf, rf / fs, oth, rf / fs - oth, tot])
    e0 = est(np.arange(n))
    # first-stage partial F
    Xf = np.c_[xi, C]; b = _ols(mi, Xf); r = mi - Xf @ b; s2 = r @ r / (n - Xf.shape[1])
    se = np.sqrt(s2 * np.linalg.inv(Xf.T @ Xf)[0, 0]); F = (b[0] / se) ** 2
    rng = np.random.default_rng(seed); B = np.array([est(rng.integers(0, n, n)) for _ in range(nboot)])
    sd = np.nanstd(B, 0)
    # use bootstrap percentile for IV CI; SE via robust (IQR/1.349) to limit weak-IV tails
    rsd = (np.nanpercentile(B, 75, 0) - np.nanpercentile(B, 25, 0)) / 1.349
    lab = ['fs', 'rf', 'iv', 'other', 'diff', 'total']
    out = dict(n=n, F=F)
    for i, l in enumerate(lab):
        out[l] = e0[i]; out[l + '_se'] = rsd[i]; out[l + '_sdboot'] = sd[i]
    # Anderson-Rubin-ish check: RF significance (reduced form)
    out['diff_p_boot'] = 2 * min(np.mean(B[:, 4] > 0), np.mean(B[:, 4] < 0))
    return out


def meta_re(b, se):
    b, se = np.asarray(b, float), np.asarray(se, float)
    ok = np.isfinite(b) & np.isfinite(se) & (se > 0); b, se = b[ok], se[ok]; k = len(b)
    if k == 0: return dict(est=np.nan, se=np.nan, lo=np.nan, hi=np.nan, p=np.nan, I2=np.nan, k=0)
    w = 1 / se ** 2; fe = (w * b).sum() / w.sum(); Q = (w * (b - fe) ** 2).sum()
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (se ** 2 + tau2); e = (wr * b).sum() / wr.sum(); s = np.sqrt(1 / wr.sum())
    return dict(est=e, se=s, lo=e - 1.96 * s, hi=e + 1.96 * s, p=2 * stats.norm.sf(abs(e / s)), I2=max(0, (Q - k + 1) / Q) if Q > 0 else 0, k=k)


def meta_table(R, keys):
    rows = []
    for k, g in R.groupby(keys):
        row = dict(zip(keys if isinstance(keys, list) else [keys], k if isinstance(k, tuple) else [k]))
        for l in ['fs', 'rf', 'iv', 'other', 'diff', 'total']:
            m = meta_re(g[l], g[l + '_se']); row[l] = m['est']; row[l + '_lo'] = m['lo']; row[l + '_hi'] = m['hi']
            if l == 'diff': row['diff_p'] = m['p']
        row['k'] = len(g); row['minF'] = g.F.min(); row['medF'] = g.F.median(); row['n'] = g.n.sum()
        rows.append(row)
    return pd.DataFrame(rows)
