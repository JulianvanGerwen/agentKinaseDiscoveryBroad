"""(a) site-level cross-study same-cell vs diff-cell replication for SAME-class pairs (inh) — is reverse a drug x cell interaction?
(b) Hijazi coverage control: reverse/primary dpct profiles recomputed on sites measured in >=90% of experiments of all 3 cells; decomposition again.
(c) Hijazi: site-level within-study, same drug diff cell vs diff class same cell (cell-vs-drug for reverse sites)."""
import numpy as np, pandas as pd
from common_i22 import get
CA = '/workspace/cache/b_I22/'; OUT = '/workspace/results/b_I22/tables/'
MP = pd.read_parquet('/workspace/results/b_I05/cache/motif_pct.parquet')
A, E = get('inh'); X = A.values.astype(np.float32); fin = ~np.isnan(X)
st = E.study.to_numpy(str); cl = E.cls.to_numpy(str); ce = E.cell.to_numpy(str); ag = E.agent.to_numpy(str)
def lor(xi, xj, sg, td=1.0, tr=0.5):
    d = sg * xi > td; r = sg * xj > tr
    if d.sum() < 10: return np.nan
    a = (d & r).sum() + .5; b = (d & ~r).sum() + .5; c = (~d & r).sum() + .5; e = (~d & ~r).sum() + .5
    return np.log(a * e / (b * c))
# (a)
I, J = np.where((st[:, None] != st[None, :]) & (cl[:, None] == cl[None, :]))
sp = np.where(st[I] < st[J], pd.Series(st[I]) + '~' + pd.Series(st[J]), pd.Series(st[J]) + '~' + pd.Series(st[I]))
df = pd.DataFrame(dict(i=I, j=J, sp=sp)); df = df[df.sp.isin(set(df.sp[ce[df.i.values] == ce[df.j.values]]))].reset_index(drop=True)
res = []
for i, j in zip(df.i, df.j):
    b = fin[:, i] & fin[:, j]
    res.append((lor(X[b, i], X[b, j], -1), lor(X[b, i], X[b, j], 1)) if b.sum() >= 500 else (np.nan, np.nan))
df[['lor_primary', 'lor_reverse']] = np.array(res); df['sc'] = ce[df.i] == ce[df.j]
g = df.groupby(['sp', 'sc'])[['lor_primary', 'lor_reverse']].mean().unstack().dropna()
out = [dict(test='site-level cross-study SAME class', dir=d_, same=g[(c, True)].mean(), diff=g[(c, False)].mean(), delta=(g[(c, True)] - g[(c, False)]).mean(),
            n_studypairs=len(g), frac_pos=((g[(c, True)] - g[(c, False)]) > 0).mean()) for d_, c in [('primary', 'lor_primary'), ('reverse', 'lor_reverse')]]
print(g.round(2).to_string())
# (c) Hijazi site-level within study
h = np.where(st == '31959955')[0]; rows = []
for a_ in h:
    for b_ in h:
        if a_ == b_: continue
        typ = 'same drug, diff cell' if (ag[a_] == ag[b_] and ce[a_] != ce[b_]) else 'diff class, same cell' if (cl[a_] != cl[b_] and ce[a_] == ce[b_]) else 'diff class, diff cell' if (cl[a_] != cl[b_]) else None
        if typ is None: continue
        if typ == 'diff class, diff cell' and np.random.default_rng(a_ * 1000 + b_).random() > 0.1: continue  # subsample
        bb = fin[:, a_] & fin[:, b_]
        rows.append(dict(type=typ, primary=lor(X[bb, a_], X[bb, b_], -1), reverse=lor(X[bb, a_], X[bb, b_], 1)))
Hs = pd.DataFrame(rows).groupby('type')[['primary', 'reverse']].agg(['mean', 'count'])
print(Hs.round(3).to_string()); Hs.to_csv(OUT + 'T4b_hijazi_sitelevel.tsv', sep='\t')
# (b) Hijazi coverage-matched profiles
cells = ['MCF7', 'HL60', 'NTERA2']
keep = np.ones(X.shape[0], bool)
for c in cells:
    keep &= fin[:, (st == '31959955') & (ce == c)].mean(1) >= 0.9
print('common sites', keep.sum())
kst = [k for k in pd.read_parquet(CA + 'profiles.parquet').columns]
Pm = MP.reindex(A.index[keep])[kst].values.astype(np.float32)
prof = {}
for cen in ['raw', 'cen']:
    Xh = X[keep][:, h].copy()
    if cen == 'cen':
        for c in cells:
            s = ce[h] == c; Xh[:, s] -= np.nanmedian(Xh[:, s], 1, keepdims=True)
    for dr, sg in [('primary', -1), ('reverse', 1)]:
        V = []; ok = []
        for k in range(len(h)):
            x = Xh[:, k]; f = ~np.isnan(x); xx = np.nan_to_num(x); bg = f & (np.abs(xx) < .3); r = f & (sg * xx > .58)
            if r.sum() < 30: continue
            V.append(Pm[r].mean(0) - Pm[bg].mean(0)); ok.append(h[k])
        Z = np.array(V); Z -= Z.mean(1, keepdims=True); Z -= Z.mean(0); Z /= np.linalg.norm(Z, axis=1, keepdims=True); R = Z @ Z.T; ok = np.array(ok)
        I2, J2 = np.triu_indices(len(ok), 1)
        sd = ag[ok][I2] == ag[ok][J2]; sc = ce[ok][I2] == ce[ok][J2]; scl = cl[ok][I2] == cl[ok][J2]
        for lab, s in [('same drug, diff cell', sd & ~sc), ('diff class, same cell', ~scl & sc), ('diff class, diff cell', ~scl & ~sc)]:
            out.append(dict(test=f'Hijazi common-site motif profile ({cen})', dir=dr, type=lab, mean_r=R[I2, J2][s].mean(), n=int(s.sum())))
T = pd.DataFrame(out); T.to_csv(OUT + 'T4_H1_controls.tsv', sep='\t', index=False); pd.set_option('display.width', 200); print(T.round(3).to_string())
