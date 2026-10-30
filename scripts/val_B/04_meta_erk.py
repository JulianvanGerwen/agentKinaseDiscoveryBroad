"""val_B step 4: (a) Fisher-z pooled rho of measured time vs scores (power); (b) is the footprint ISCH score 'another MAPK readout'? coherence after regressing an ERK activation-loop readout (excl. T202, which is in ISCH set) vs ISCH score; partial correlations."""
import sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/audit_I10'); sys.path.insert(0, '/workspace/scripts/val_B')
from common import *
from scipy import stats
RV = '/workspace/results/val_B/'
C = pd.read_csv(RV + 'T2_correlations.tsv', sep='\t'); m = C[C.test == 'measured time']
rows = []
for v, d in m.groupby('var'):
    z = np.arctanh(d.rho.values); w = d.n.values - 3; zbar = (z * w).sum() / w.sum(); se = 1 / np.sqrt(w.sum())
    Q = (w * (z - zbar) ** 2).sum(); k = len(z)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); wr = 1 / (1 / w + tau2); zr = (z * wr).sum() / wr.sum(); ser = 1 / np.sqrt(wr.sum())
    rows.append(dict(var=v, k=k, N=d.n.sum(), fixed_rho=np.tanh(zbar), lo=np.tanh(zbar - 1.96 * se), hi=np.tanh(zbar + 1.96 * se), p=2 * stats.norm.sf(abs(zbar / se)),
                     re_rho=np.tanh(zr), re_lo=np.tanh(zr - 1.96 * ser), re_hi=np.tanh(zr + 1.96 * ser), I2=max(0, (Q - (k - 1)) / Q) if Q > 0 else 0))
M = pd.DataFrame(rows); M.to_csv(RV + 'T4_meta_time_corr.tsv', sep='\t', index=False); print(M.round(3).to_string())
SETS = {'MEK': 'KI_MEK', 'EGF': 'EGF', 'TNF': 'TNF', 'DNAdmg': 'DNA_damage'}
ISCH = ischaemia_set(200)
ERK = ['P28482_T185', 'P28482_Y187', 'P27361_Y204']          # T202 omitted (member of ISCH set)
def resid(Z, v):
    v = (v - v.mean()) / v.std(); v = np.nan_to_num(v.values); X = Z.values; M_ = np.isfinite(X)
    with np.errstate(all='ignore'): b = np.nansum(np.where(M_, X, 0) * v, 1) / (M_ * v ** 2).sum(1)
    return standardise(pd.DataFrame(np.where(M_, X - b[:, None] * v, np.nan), index=Z.index, columns=Z.columns))
rows = []
for coh in COH:
    Z = standardise(cohort(coh)); xr = cohort(coh, mincov=0.0); xr = xr.loc[[s for s in ERK if s in xr.index]]
    erk = standardise(xr).mean(0).reindex(Z.columns); isc = score(Z, ISCH[ISCH.index.isin(Z.index)])
    sc = {k: score(Z, fp_set(s)[lambda q: q.index.isin(Z.index)]) for k, s in SETS.items()}
    # drop footprint-site overlap with the ERK sites (none expected) -> regress
    for lab, cov in [('none', None), ('footprint ISCH score', isc), ('ERK activation-loop readout (no T202)', erk)]:
        Zx = Z if cov is None else resid(Z, cov); Cx = Coh(Zx)
        for k, s in SETS.items():
            sg = fp_set(s); sg = sg[sg.index.isin(Cx.idx)]; ii = np.array([Cx.idx[q] for q in sg.index])
            rows.append(dict(cohort=coh, set=k, adj=lab, coherence=Cx.coherence(ii, sg.values)[0]))
    # partial Spearman MEK~ISCH given ERK
    for k in ['MEK', 'EGF', 'TNF']:
        d = pd.DataFrame({'a': sc[k], 'b': isc, 'e': erk}).dropna(); rk = d.rank()
        res = lambda y: y - np.polyval(np.polyfit(rk.e, y, 1), rk.e)
        rows.append(dict(cohort=coh, set=k, adj='partial rho(score, ISCH | ERK readout)', coherence=stats.pearsonr(res(rk.a), res(rk.b))[0]))
        rows.append(dict(cohort=coh, set=k, adj='rho(score, ERK readout)', coherence=stats.spearmanr(d.a, d.e)[0]))
    rows.append(dict(cohort=coh, set='ISCH', adj='rho(ISCH, ERK readout)', coherence=stats.spearmanr(*pd.concat([isc, erk], axis=1).dropna().values.T)[0]))
    print(coh, flush=True)
R = pd.DataFrame(rows); R.to_csv(RV + 'T5_erk_vs_isch_adjust.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(R.groupby(['adj', 'set']).coherence.median().unstack(1).round(3).to_string())
