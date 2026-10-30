"""val_B: measured (GDC) cold-ischaemia time vs footprint scores, and coherence after adjustment. Reuses audit_I10/common.py (cohort, standardise, Coh, score, fp_set, ischaemia_set)."""
import sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/audit_I10')
from common import *
from scipy import stats
RV = '/workspace/results/val_B/'
g = pd.read_csv(RV + 'gdc_ischaemia_all.tsv', sep='\t', index_col=0)
g.loc[g.excision_min > 1440, 'excision_min'] = np.nan          # 1 implausible value (525618 min) in LSCC
VE = pd.read_csv('/workspace/results/val_E/ischaemia_scores.tsv', sep='\t', index_col=0); VE = VE[VE.variant == 'input_tc_padj']
ISCH = ischaemia_set(200)
SETS = {'MEK': 'KI_MEK', 'EGF': 'EGF', 'TNF': 'TNF', 'DNAdmg': 'DNA_damage'}
ERK = ['P28482_T185', 'P28482_Y187', 'P27361_T202', 'P27361_Y204']
rng = np.random.default_rng(5)

def resid(Z, covs):
    """regress covariate matrix (samples x k, no NaN) out of every site (observed values), re-standardise."""
    C = np.column_stack([np.ones(len(covs))] + [covs[:, j] for j in range(covs.shape[1])])
    X = Z.values; out = np.full_like(X, np.nan); M = np.isfinite(X)
    for i in range(X.shape[0]):
        m = M[i]
        if m.sum() < C.shape[1] + 5: continue
        b = np.linalg.lstsq(C[m], X[i, m], rcond=None)[0]; out[i, m] = X[i, m] - C[m] @ b
    return standardise(pd.DataFrame(out, index=Z.index, columns=Z.columns))

def zs(v): v = np.asarray(v, float); return (v - v.mean()) / v.std()
corr_rows, coh_rows = [], []
for coh in COH:
    x = cohort(coh); Z = standardise(x); tum = list(Z.columns); cid = [t.rsplit('_', 1)[0] for t in tum]
    t = pd.Series(g.excision_min.reindex(cid).values, index=tum)
    have = t.notna().sum()
    sc = {k: score(Z, fp_set(s)[lambda q: q.index.isin(Z.index)]) for k, s in SETS.items()}
    sc['ISCH_audit'] = score(Z, ISCH[ISCH.index.isin(Z.index)])
    sc['ISCH_valE'] = VE.isch.reindex(tum)
    sc['ISCH_valE_noPath'] = VE.isch_noPath.reindex(tum)
    xr = cohort(coh, mincov=0.0); xr = xr.loc[[s for s in ERK if s in xr.index]]
    zz = standardise(xr); sc['ERK_loop'] = zz.mean(0).reindex(tum) if len(zz) else pd.Series(np.nan, index=tum)
    S = pd.DataFrame(sc).reindex(tum)
    # ---- task 2: correlations
    if have >= 20:
        for k in S:
            d = pd.concat([S[k], t], axis=1).dropna()
            rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1]); rl = stats.pearsonr(d.iloc[:, 0], np.log(d.iloc[:, 1]))[0]
            corr_rows.append(dict(cohort=coh, var=k, test='measured time', n=len(d), rho=rho, p=p, r_logt=rl))
        for k in ['MEK', 'EGF', 'TNF', 'DNAdmg', 'ERK_loop']:
            for ik in ['ISCH_audit']:
                d = S[[k, ik]].dropna(); rho, p = stats.spearmanr(d[k], d[ik]); corr_rows.append(dict(cohort=coh, var=k, test=f'vs {ik} (all tumours)', n=len(d), rho=rho, p=p))
                d = S[[k, ik]].join(t.rename('t')).dropna(); rho, p = stats.spearmanr(d[k], d[ik])
                corr_rows.append(dict(cohort=coh, var=k, test=f'vs {ik} (time-known subset)', n=len(d), rho=rho, p=p))
                # partial Spearman of footprint vs ISCH given log t
                rk = d.rank(); rk['lgt'] = np.log(d.t)
                res = lambda y: y - np.polyval(np.polyfit(rk.lgt, y, 1), rk.lgt)
                pr = stats.pearsonr(res(rk[k]), res(rk[ik]))
                corr_rows.append(dict(cohort=coh, var=k, test=f'partial vs {ik} | log t', n=len(d), rho=pr[0], p=pr[1]))
        # variance of ISCH score explained by time
        d = S[['ISCH_audit']].join(t.rename('t')).dropna(); T = np.column_stack([d.t, np.log(d.t)])
        r2 = 1 - np.var(d.ISCH_audit - np.column_stack([np.ones(len(d)), T]) @ np.linalg.lstsq(np.column_stack([np.ones(len(d)), T]), d.ISCH_audit, rcond=None)[0]) / np.var(d.ISCH_audit)
        corr_rows.append(dict(cohort=coh, var='ISCH_audit', test='R2 of score on t+log t', n=len(d), rho=r2))
    # ---- task 3: coherence
    def coh_of(Zx, label, n):
        C = Coh(Zx)
        for k, s in SETS.items():
            sg = fp_set(s); sg = sg[sg.index.isin(C.idx)]; ii = np.array([C.idx[q] for q in sg.index])
            coh_rows.append(dict(cohort=coh, set=k, adj=label, n_tumours=n, coherence=C.coherence(ii, sg.values)[0]))
    if have >= 20:
        keep = t.notna().values; Zs = standardise(Z.loc[:, keep]); ts = t[keep].values; n = keep.sum()
        coh_of = coh_of
        coh_of(Zs, 'none', n)
        coh_of(resid(Zs, np.column_stack([ts, np.log(ts)])), 'measured t + log t', n)
        coh_of(resid(Zs, np.column_stack([np.log(ts)])), 'measured log t', n)
        coh_of(resid(Zs, zs(S.ISCH_audit[keep]).reshape(-1, 1)) if S.ISCH_audit[keep].notna().all() else Zs, 'footprint ISCH score', n)
        vv = S.ISCH_valE_noPath[keep]
        if vv.notna().all(): coh_of(resid(Zs, zs(vv).reshape(-1, 1)), 'val_E ISCH noPath score', n)
        coh_of(resid(Zs, np.column_stack([zs(S.ISCH_audit[keep]), ts, np.log(ts)])) if S.ISCH_audit[keep].notna().all() else Zs, 'ISCH score + t + log t', n)
        for b in range(3):
            pt = rng.permutation(ts); coh_of(resid(Zs, np.column_stack([pt, np.log(pt)])), f'PERMUTED t + log t #{b}', n)
        # random 2-df covariate control matched to footprint score (permuted ISCH score + its square)
    # full cohort (all tumours): none vs ISCH score (audit reproduction), incl. cohorts without times
    coh_of(Z, 'ALL-tumours none', len(tum))
    if S.ISCH_audit.notna().all(): coh_of(resid(Z, zs(S.ISCH_audit).reshape(-1, 1)), 'ALL-tumours footprint ISCH score', len(tum))
    print(coh, have, flush=True)
C = pd.DataFrame(corr_rows); C.to_csv(RV + 'T2_correlations.tsv', sep='\t', index=False)
H = pd.DataFrame(coh_rows); H.to_csv(RV + 'T3_coherence.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
m = C[C.test == 'measured time']
print(m.pivot(index='var', columns='cohort', values='rho').round(2).to_string())
print(m.groupby('var').rho.agg(['median', 'min', 'max']).round(2))
for tst in ['vs ISCH_audit (all tumours)', 'vs ISCH_audit (time-known subset)', 'partial vs ISCH_audit | log t']:
    print(tst); print(C[C.test == tst].pivot(index='var', columns='cohort', values='rho').round(2).to_string())
print(C[C.test.str.startswith('R2')][['cohort', 'rho']].round(3).to_string())
print(H.pivot_table(index=['set', 'adj'], columns='cohort', values='coherence').round(3).to_string())
print(H.groupby(['set', 'adj']).coherence.median().unstack(0).round(3))
