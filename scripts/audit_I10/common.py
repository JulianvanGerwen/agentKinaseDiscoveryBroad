"""audit_I10 shared helpers (independent re-implementation; reads only infra assets + b_I10 footprint table T1)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
sys.path.insert(0, '/workspace/scripts/infra_compendium'); import load
CA = '/workspace/cache/audit_I10/'; RA = '/workspace/results/audit_I10/'
COH = cl.COHORTS
FP = pd.read_csv('/workspace/results/b_I10/tables/T1_footprints.tsv', sep='\t')

def fp_set(g):
    d = FP[FP.group == g]; return pd.Series(np.sign(d.w.values), index=d.site.values)

_S = None
def site_annot():
    global _S
    if _S is None:
        S = load.sites().set_index('site')
        S['pro1'] = S.sequence_window.str[8] == 'P'
        _S = S
    return _S

def psp_sets(minn=15):
    S = site_annot(); ks = {}
    for s, v in S.psp_kinases_human.dropna().items():
        for k in v.split(';'): ks.setdefault(k, []).append(s)
    return {k: pd.Series(1.0, index=v) for k, v in ks.items() if len(v) >= minn}

def cohort(coh, variant='input_tc_padj', mincov=0.25):
    x = cl.matrix(coh, variant); tum = [s for s in cl.tumour_ids(coh) if s in x.columns]
    x = x[tum]; x = x[(x.notna().mean(1) >= mincov) & (x.notna().sum(1) >= 10)]
    return x.astype(np.float64)

def standardise(x):
    mu = x.mean(1); sd = x.std(1)
    return x.sub(mu, axis=0).div(sd.replace(0, np.nan), axis=0)

class Coh:
    """Per-cohort standardised matrix with fast pairwise-complete correlations (sites standardised over observed values)."""
    def __init__(self, Z):
        self.Z = Z; self.idx = {s: i for i, s in enumerate(Z.index)}
        self.M = np.isfinite(Z.values).astype(np.float64); self.X0 = np.nan_to_num(Z.values)
        self.prot = np.array([s.split('_')[0] for s in Z.index])
    def corr(self, ii):
        num = self.X0[ii] @ self.X0[ii].T; n = self.M[ii] @ self.M[ii].T
        with np.errstate(all='ignore'): r = num / n
        r[n < 20] = np.nan; return r
    def coherence(self, ii, sg, extra_mask=None):
        r = self.corr(ii); P = np.outer(sg, sg) * r
        m = np.triu(np.ones_like(r, bool), 1) & (self.prot[ii][:, None] != self.prot[ii][None, :]) & np.isfinite(P)
        if extra_mask is not None: m &= extra_mask
        return P[m].mean() if m.sum() >= 10 else np.nan, int(m.sum())

def residualise_pcs(Z, k, mincov_pca=0.5):
    """Remove top-k sample PCs (computed on sites with >=50% coverage, mean-imputed) from every site (fit on observed values)."""
    if k == 0: return Z
    A = Z[Z.notna().mean(1) >= mincov_pca].fillna(0).values
    U, s, Vt = np.linalg.svd(A - A.mean(1, keepdims=True), full_matrices=False)
    F = Vt[:k].T                                     # samples x k
    X = Z.values.copy(); out = np.full_like(X, np.nan)
    M = np.isfinite(X)
    # group sites by missingness pattern is overkill; loop with lstsq on observed samples
    for i in range(X.shape[0]):
        m = M[i]
        if m.sum() < k + 5: continue
        Fm = np.column_stack([np.ones(m.sum()), F[m]])
        b = np.linalg.lstsq(Fm, X[i, m], rcond=None)[0]
        out[i, m] = X[i, m] - Fm @ b
    R = pd.DataFrame(out, index=Z.index, columns=Z.columns)
    return standardise(R), F, s

def score(Z, sg, minsites=10):
    """Footprint score per tumour: signed mean of standardised site values over observed footprint sites."""
    sg = sg[sg.index.isin(Z.index)]
    A = Z.loc[sg.index].mul(sg, axis=0)
    n = A.notna().sum(0)
    return A.mean(0).where(n >= minsites)

def ischaemia_set(k=200):
    """Top-k sites of the Mertins 2014 cold-ischaemia experiments in the compendium (sign-consistent, >=2 measured)."""
    W = pd.read_parquet('/workspace/cache/b_I10/fp_sites.parquet')
    I = W[[c for c in W.columns if c.startswith('24719451|') and 'ischaemia' in c]]
    n = I.notna().sum(1); mu = I.mean(1); cons = (np.sign(I).eq(np.sign(mu), axis=0) | I.isna()).all(1)
    r = mu[(n >= 2) & cons].abs().sort_values(ascending=False).index[:k]
    return pd.Series(np.sign(mu[r]).values, index=r)
