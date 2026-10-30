"""b_Q8 shared helpers: tumour vs NAT (normal adjacent tissue) preparation with symmetric, within-group processing.
Per cohort: site universe = sites with >=25% coverage and >=10 observations in BOTH tumours and NATs (raw log2 phospho).
Within each group separately: (optional) per-site OLS on own protein (phospho ~ a + b*protein; residual), (optional) extra covariates
(composition / ischaemia), (optional) top-k sample PCs removed; then site-wise z over observed values."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
sys.path.insert(0, '/workspace/scripts/infra_compendium'); import load
RQ = '/workspace/results/b_Q8/'; TQ = RQ + 'tables/'; FQ = RQ + 'figures/'; CQ = '/workspace/cache/b_Q8/'
NATC = ['CCRCC', 'COAD', 'HNSCC', 'LSCC', 'LUAD', 'PDAC', 'OV', 'UCEC']      # cohorts with NATs
MAINC = ['CCRCC', 'COAD', 'HNSCC', 'LSCC', 'LUAD', 'PDAC']                   # NAT n >= 44
FP = pd.read_csv('/workspace/results/b_I10/tables/T1_footprints.tsv', sep='\t')
GROUPS = sorted(FP.group.unique())

def fp_set(g):
    d = FP[FP.group == g]; return pd.Series(np.sign(d.w.values), index=d.site.values)

_S = None
def site_annot():
    global _S
    if _S is None:
        S = load.sites().set_index('site'); S['pro1'] = S.sequence_window.str[8] == 'P'; _S = S
    return _S

def psp_sets(minn=15):
    S = site_annot(); ks = {}
    for s, v in S.psp_kinases_human.dropna().items():
        for k in v.split(';'): ks.setdefault(k, []).append(s)
    return {k: sorted(set(v)) for k, v in ks.items() if len(set(v)) >= minn}

def groups(coh):
    s = cl.samples(); s = s[(s.cohort == coh) & s.has_phospho]
    T = s[s.sample_type == 'Tumor']; N = s[s.sample_type == 'Normal']
    paired = sorted(set(T.case_id) & set(N.case_id))
    return list(T.sample_id), list(N.sample_id), paired

def load_raw(coh):
    ph = cl.matrix(coh, 'phospho').astype(np.float64); pr = cl.matrix(coh, 'protein').astype(np.float64)
    T, N, paired = groups(coh)
    T = [c for c in T if c in ph.columns]; N = [c for c in N if c in ph.columns]
    ok = lambda X: (X.notna().mean(1) >= 0.25) & (X.notna().sum(1) >= 10)
    keep = ok(ph[T]) & ok(ph[N])
    ph = ph[keep]
    return ph, pr, T, N, paired

def standardise(x):
    return x.sub(x.mean(1), axis=0).div(x.std(1).replace(0, np.nan), axis=0)

def _resid(X, F, minobs):
    """Row-wise OLS of X (sites x samples, NaN allowed) on sample covariates F (samples x p, no NaN); returns residuals."""
    out = np.full_like(X, np.nan); M = np.isfinite(X)
    pats, inv = np.unique(M, axis=0, return_inverse=True)
    for k in range(len(pats)):
        m = pats[k]
        if m.sum() < F.shape[1] + minobs: continue
        rows = np.where(inv.ravel() == k)[0]
        Fm = np.column_stack([np.ones(m.sum()), F[m]])
        B = np.linalg.lstsq(Fm, X[np.ix_(rows, np.where(m)[0])].T, rcond=None)[0]
        out[np.ix_(rows, np.where(m)[0])] = X[np.ix_(rows, np.where(m)[0])] - (Fm @ B).T
    return out

def protein_adjust(ph, pr, minobs=10):
    """Per-site OLS phospho ~ protein within the given samples; sites without protein / too few pairs -> dropped."""
    u = np.array([s.split('_')[0] for s in ph.index]); P = pr.reindex(index=u, columns=ph.columns).values
    X = ph.values; ok = np.isfinite(X) & np.isfinite(P); n = ok.sum(1)
    Xm = np.where(ok, X, np.nan); Pm = np.where(ok, P, np.nan)
    with np.errstate(all='ignore'):
        mx = np.nanmean(Xm, 1, keepdims=True); mp = np.nanmean(Pm, 1, keepdims=True)
        b = np.nansum((Xm - mx) * (Pm - mp), 1) / np.nansum((Pm - mp) ** 2, 1)
    R = Xm - mx - b[:, None] * (Pm - mp)
    R[n < minobs] = np.nan
    return pd.DataFrame(R, index=ph.index, columns=ph.columns), b

def global_pcs(Z, k, mincov=0.5):
    A = Z[Z.notna().mean(1) >= mincov].fillna(0).values
    U, s, Vt = np.linalg.svd(A - A.mean(1, keepdims=True), full_matrices=False)
    return Vt[:k].T, s ** 2 / (s ** 2).sum()

def process(ph, pr, samples, padj=True, k=0, cov=None, minobs=10):
    """Returns standardised site x sample matrix for one group. cov: DataFrame samples x p (no NaN) regressed out."""
    X = ph[samples]
    if padj: X, _ = protein_adjust(X, pr, minobs)
    X = X[X.notna().sum(1) >= minobs]
    if cov is not None:
        F = cov.loc[samples].values; F = (F - F.mean(0)) / F.std(0)
        X = pd.DataFrame(_resid(X.values, F, minobs), index=X.index, columns=X.columns)
    Z = standardise(X)
    if k > 0:
        F, _ = global_pcs(Z, k)
        Z = standardise(pd.DataFrame(_resid(Z.values, F, minobs), index=Z.index, columns=Z.columns))
    return Z.dropna(how='all')

class Coh:
    """Fast pairwise-complete correlations on a standardised matrix (observed values standardised)."""
    def __init__(self, Z, minpair=15):
        self.Z = Z; self.idx = {s: i for i, s in enumerate(Z.index)}; self.minpair = minpair
        self.M = np.isfinite(Z.values).astype(np.float64); self.X0 = np.nan_to_num(Z.values)
        self.prot = np.array([s.split('_')[0] for s in Z.index])
    def corr(self, ii, jj=None):
        jj = ii if jj is None else jj
        num = self.X0[ii] @ self.X0[jj].T; n = self.M[ii] @ self.M[jj].T
        with np.errstate(all='ignore'): r = num / (n - 1)
        r = np.clip(r, -1, 1); r[n < self.minpair] = np.nan; return r
    def coherence(self, ii, sg):
        r = self.corr(ii); P = np.outer(sg, sg) * r
        m = np.triu(np.ones_like(r, bool), 1) & (self.prot[ii][:, None] != self.prot[ii][None, :]) & np.isfinite(P)
        return (P[m].mean() if m.sum() >= 10 else np.nan), int(m.sum())

def score(Z, sg, minsites=5):
    sg = sg[sg.index.isin(Z.index)]
    A = Z.loc[sg.index].mul(sg, axis=0); n = A.notna().sum(0)
    return A.mean(0).where(n >= minsites)

def composition_markers():
    return pd.read_csv(TQ + 'T0_composition_markers.tsv', sep='\t')

def ischaemia_set():
    return pd.read_csv(CQ + 'ischaemia_fp.tsv', sep='\t', index_col=0).iloc[:, 0]
