"""val_E step 3: KRAS offset (b_I24 T4 model, RTK>ERK, tc_padj, mut) and LUAD EGFR-mut offset, +/- ischaemia covariate."""
import sys, numpy as np, pandas as pd, warnings, statsmodels.api as sm; warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/audit_I24'); from common_a import cl, zs, meta_re, COH
OUT = '/workspace/results/val_E/'
ST = pd.read_pickle('/workspace/cache/val_E/st_padj.pkl')
def fit(y, x, l, C, extra_cols=0):
    """exact b_I24 04_modifiers 'main' model: y~x+L+x:L+C+x:C[:2]; Huber RLM; extra C columns enter additively."""
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(l) & np.all(np.isfinite(C), 1)
    if ok.sum() < 30 or l[ok].sum() < 5 or (1 - l[ok]).sum() < 5: return None
    y_, x_, l_, C_ = zs(y[ok]), zs(x[ok]), l[ok], np.apply_along_axis(zs, 0, C[ok])
    cols = [x_, l_, x_ * l_] + [C_[:, j] for j in range(C_.shape[1])] + [x_ * C_[:, j] for j in range(2)]
    r = sm.RLM(y_, sm.add_constant(np.column_stack(cols)), M=sm.robust.norms.HuberT()).fit()
    return dict(n=int(ok.sum()), n_mut=int(l[ok].sum()), b_lev=r.params[2], se_lev=r.bse[2], b_int=r.params[3], se_int=r.bse[3])
def mutvec(c, D, genes):
    mb = cl.matrix(c, 'mut_binary'); v = np.zeros(len(D))
    for g in genes:
        if g in mb.index: v = np.maximum(v, mb.loc[g].reindex(D.case_id).fillna(0).values)
    v = v.astype(float); v[~D.case_id.isin(mb.columns).values] = np.nan; return v
rows = []
for c in COH:
    D = ST[('tc_padj', c)]
    if D.RTK.notna().sum() < 30: continue
    k = mutvec(c, D, ['KRAS'])
    for nm, cv in [('base', ['purity', 'prolif']), ('+isch', ['purity', 'prolif', 'isch'])]:
        r = fit(D.ERK.values.astype(float), D.RTK.values.astype(float), k, D[cv].values.astype(float))
        if r: rows.append(dict(test='KRASmut RTK>ERK offset', cohort=c, spec=nm, **r))
    if c == 'LUAD':
        e = mutvec(c, D, ['EGFR'])
        for nm, cv in [('base', ['purity', 'prolif']), ('+isch', ['purity', 'prolif', 'isch'])]:
            r = fit(D.ERK.values.astype(float), D.RTK.values.astype(float), e, D[cv].values.astype(float))
            if r: rows.append(dict(test='EGFRmut RTK>ERK offset (LUAD, vs all non-mut)', cohort=c, spec=nm, **r))
# LUAD: EGFR-mut vs driver-negative, ERK pTpY level (Huber, standardised)
D = ST[('tc_padj', 'LUAD')]; mb = cl.matrix('LUAD', 'mut_binary')
drv = mutvec('LUAD', D, ['KRAS', 'EGFR', 'BRAF', 'NRAS', 'HRAS', 'ERBB2', 'MET', 'RIT1', 'MAP2K1']); egfr = mutvec('LUAD', D, ['EGFR'])
oth = mutvec('LUAD', D, ['KRAS', 'BRAF', 'NRAS', 'HRAS', 'ERBB2', 'MET', 'RIT1', 'MAP2K1'])
grp = np.where(egfr == 1, 1., np.where(drv == 0, 0., np.nan)); grp[(egfr == 1) & (oth == 1)] = np.nan
print('LUAD EGFR-mut n', int((grp == 1).sum()), 'driver-neg n', int((grp == 0).sum()))
for nm, cv in [('base', ['purity', 'prolif']), ('+isch', ['purity', 'prolif', 'isch'])]:
    d = D.assign(g=grp)[['ERK', 'g'] + cv].dropna(); y = zs(d.ERK.values); C = np.apply_along_axis(zs, 0, d[cv].values.astype(float))
    r = sm.RLM(y, sm.add_constant(np.c_[d.g.values, C]), M=sm.robust.norms.HuberT()).fit()
    rows.append(dict(test='LUAD ERK level: EGFR-mut vs driver-neg', cohort='LUAD', spec=nm, n=len(d), n_mut=int(d.g.sum()), b_lev=r.params[1], se_lev=r.bse[1]))
R = pd.DataFrame(rows); R.to_csv(OUT + 'offsets_by_cohort.tsv', sep='\t', index=False)
K = R[R.test.str.startswith('KRAS')]; out = []
for s, g in K.groupby('spec', sort=False):
    m = meta_re(g.b_lev, g.se_lev); out.append(dict(test='KRASmut RTK>ERK offset (meta)', spec=s, k=len(g), est=m['est'], lo=m['lo'], hi=m['hi'], p=m['p'], n_mut=g.n_mut.sum()))
pd.DataFrame(out).to_csv(OUT + 'offsets_meta.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(pd.DataFrame(out).round(3).to_string()); print(R.round(3).to_string())
# isch vs KRAS within LUAD/COAD/PDAC/UCEC: difference in isch score (SD units) KRAS-mut vs wt
for c in COH:
    D = ST[('tc_padj', c)]; k = mutvec(c, D, ['KRAS']); d = pd.DataFrame({'i': zs(D.isch.values), 'k': k}).dropna()
    if d.k.sum() >= 5: 
        a, b = d.i[d.k == 1], d.i[d.k == 0]; print(c, 'KRAS-mut isch diff (SD) %.2f  n_mut %d  p=%.2g' % (a.mean() - b.mean(), len(a), __import__('scipy.stats').stats.mannwhitneyu(a, b).pvalue))
