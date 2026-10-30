"""Exploratory genome-wide scan of CNV/RNA gain modifiers (OLS via Frisch-Waugh, vectorized), RE meta, FDR, >=3-cohort replication.
Permutation control: lesion sample labels permuted within cohort (same perm across genes). Usage: 05_scan.py [n_perm]"""
from common import *
from scipy import stats
from statsmodels.stats.multitest import multipletests
import sys
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 3
A = load_stages('tc_padj'); COV = ['purity', 'prolif']; rng = np.random.default_rng(7)

def prep_feats(c, D, t):
    M = cl.matrix(c, 'cnv_log2' if t == 'cnv' else 'rna')
    M = M.reindex(columns=(D.case_id.values if t == 'cnv' else D.index))
    M.columns = D.index
    M = M[M.notna().mean(1) >= 0.8]
    if t == 'rna':
        M = M[(M.mean(1) > M.mean(1).quantile(0.3)) & (M.std(1) > 0.1)]
    else:
        M = M[M.std(1) > 0.05]
    return M

def scan(y, x, C, F):
    """y,x (n,), C (n,k), F (n,G) -> b_int, se_int, b_lev, se_lev arrays (G,)"""
    n = len(y); y, x = zs(y), zs(x); C = np.apply_along_axis(zs, 0, C)
    F = (F - np.nanmean(F, 0)) / np.nanstd(F, 0); F = np.where(np.isfinite(F), F, 0)
    B = np.column_stack([np.ones(n), x, C, x[:, None] * C])
    Q, _ = np.linalg.qr(B)
    res = lambda Z: Z - Q @ (Q.T @ Z)
    ry = res(y[:, None])[:, 0]; L = res(F); XL = res(x[:, None] * F)
    a11 = (L * L).sum(0); a22 = (XL * XL).sum(0); a12 = (L * XL).sum(0)
    b1 = L.T @ ry; b2 = XL.T @ ry; det = a11 * a22 - a12**2
    bl = (a22 * b1 - a12 * b2) / det; bi = (a11 * b2 - a12 * b1) / det
    rss = (ry**2).sum() - (bl * b1 + bi * b2); df = n - B.shape[1] - 2; s2 = rss / df
    return bi, np.sqrt(s2 * a11 / det), bl, np.sqrt(s2 * a22 / det)

def meta_vec(B, S):
    W = 1 / S**2; W = np.where(np.isfinite(W), W, 0); B0 = np.where(np.isfinite(B), B, 0)
    k = (W > 0).sum(1); sw = W.sum(1); fe = (W * B0).sum(1) / sw
    Qs = (W * (B0 - fe[:, None])**2).sum(1); tau2 = np.maximum(0, (Qs - (k - 1)) / (sw - (W**2).sum(1) / sw))
    Wr = np.where(W > 0, 1 / (1 / np.where(W > 0, W, 1) + tau2[:, None]), 0)
    est = (Wr * B0).sum(1) / Wr.sum(1); se = np.sqrt(1 / Wr.sum(1))
    return est, se, 2 * stats.norm.sf(abs(est / se)), k

results = []; permsum = []
for t in ['rna', 'cnv']:
    for e in EDGES:
        xs, ys = e; per = {}
        for c in COHORTS:
            D = A[A.cohort == c]; D = D[D[[xs, ys] + COV].notna().all(1)]
            if len(D) < 30: continue
            F = prep_feats(c, D, t)
            out = {}
            for p in range(NPERM + 1):
                Fv = F.values.T if p == 0 else F.values.T[rng.permutation(len(D))]
                bi, si, bl, sl = scan(D[ys].values, D[xs].values, D[COV].values, Fv)
                out[p] = pd.DataFrame({'bi': bi, 'si': si, 'bl': bl, 'sl': sl}, index=F.index)
            per[c] = out
        if len(per) < 5: continue
        for p in range(NPERM + 1):
            genes = sorted(set().union(*[per[c][p].index for c in per]))
            Bi = np.column_stack([per[c][p].bi.reindex(genes) for c in per]); Si = np.column_stack([per[c][p].si.reindex(genes) for c in per])
            Bl = np.column_stack([per[c][p].bl.reindex(genes) for c in per]); Sl = np.column_stack([per[c][p].sl.reindex(genes) for c in per])
            ei, sei, pi, k = meta_vec(Bi, Si); el, sel, pl, _ = meta_vec(Bl, Sl)
            Zi = Bi / Si; Zl = Bl / Sl
            repi = ((np.sign(Zi) == np.sign(ei)[:, None]) & (np.abs(Zi) > 1.96)).sum(1)
            repl = ((np.sign(Zl) == np.sign(el)[:, None]) & (np.abs(Zl) > 1.96)).sum(1)
            df = pd.DataFrame({'type': t, 'edge': ename(e), 'gene': genes, 'k': k, 'int': ei, 'int_se': sei, 'int_p': pi, 'int_rep': repi,
                               'lev': el, 'lev_se': sel, 'lev_p': pl, 'lev_rep': repl})
            df = df[df.k >= 5]
            df['int_q'] = multipletests(df.int_p, method='fdr_bh')[1]; df['lev_q'] = multipletests(df.lev_p, method='fdr_bh')[1]
            hi = ((df.int_q < 0.05) & (df.int_rep >= 3)).sum(); hl = ((df.lev_q < 0.05) & (df.lev_rep >= 3)).sum()
            permsum.append(dict(type=t, edge=ename(e), perm=p, n_tests=len(df), int_hits=hi, lev_hits=hl,
                                int_lambda=np.median(stats.chi2.isf(df.int_p, 1)) / 0.455, lev_lambda=np.median(stats.chi2.isf(df.lev_p, 1)) / 0.455))
            if p == 0: results.append(df)
        print(permsum[-NPERM - 1:], flush=True)
R = pd.concat(results); R.to_parquet(f'{OUT}/tables/T5_scan_all.parquet')
P = pd.DataFrame(permsum); P.to_csv(f'{OUT}/tables/T5_scan_perm_summary.tsv', sep='\t', index=False)
H = R[((R.int_q < 0.05) & (R.int_rep >= 3)) | ((R.lev_q < 0.05) & (R.lev_rep >= 3))]
H.to_csv(f'{OUT}/tables/T5_scan_hits.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(P.to_string())
print(R.sort_values('int_p').head(20).to_string())
