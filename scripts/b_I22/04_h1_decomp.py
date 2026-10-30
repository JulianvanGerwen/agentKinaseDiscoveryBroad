"""H1 follow-ups. (i) Hijazi decomposition (3 cells x 58 drugs, shared DMSO per cell): profile r for same drug/diff cell, diff drug/same cell, diff/diff.
(ii) Site-level: reverse (and primary) responders in experiment i (|x|>1) replicate in j (x>0.5 in same direction) — LOR, matched within study pair,
same cell vs diff cell, cross-study + different class only; inh-inh and stim-stim pairs; block permutation of cell labels over units."""
import numpy as np, pandas as pd
from common_i22 import get
CA = '/workspace/cache/b_I22/'; OUT = '/workspace/results/b_I22/tables/'
P = pd.read_parquet(CA + 'profiles.parquet'); M = pd.read_csv(CA + 'profiles_meta.tsv', sep='\t', dtype={'study': str})
rows = []
for defn in ['dpct', 'lor90']:
    for cen in ['raw', 'cen']:
        for dr in ['primary', 'reverse']:
            ix = M.index[(M.defn == defn) & (M.cen == cen) & (M.dir == dr) & (M.kind == 'inh')].values
            m = M.loc[ix].reset_index(drop=True); Z = P.loc[ix].values.astype(float); Z -= Z.mean(1, keepdims=True); Z -= Z.mean(0); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
            h = m.index[m.study == '31959955'].values; R = Z[h] @ Z[h].T; mh = m.loc[h]
            I, J = np.triu_indices(len(h), 1)
            sd = mh.agent.values[I] == mh.agent.values[J]; sc = mh.cell.values[I] == mh.cell.values[J]; scl = mh.cls.values[I] == mh.cls.values[J]
            for lab, s in [('same drug, diff cell', sd & ~sc), ('diff class, same cell', ~scl & sc), ('diff class, diff cell', ~scl & ~sc)]:
                rows.append(dict(defn=defn, cen=cen, dir=dr, type=lab, mean_r=R[I, J][s].mean(), n=int(s.sum())))
import sys
H = pd.DataFrame(rows).pivot_table(index=['defn', 'cen', 'dir'], columns='type', values='mean_r')
H.to_csv(OUT + 'T2_hijazi_decomposition.tsv', sep='\t'); print(H.round(3).to_string())

# ---- site level
def site_level(kind, thr_d=1.0, thr_r=0.5, nperm=500):
    A, E = get(kind); X = A.values.astype(np.float32); fin = ~np.isnan(X)
    prim = -1 if kind == 'inh' else 1
    E['unit'] = E.study + '|' + E.cell
    st = E.study.to_numpy(str); cl = E.cls.to_numpy(str); ce = E.cell.to_numpy(str)
    # candidate pairs: cross-study, diff class, study pair that contains >=1 same-cell pair
    I, J = np.where((st[:, None] != st[None, :]) & (cl[:, None] != cl[None, :]))
    sp = np.where(st[I] < st[J], pd.Series(st[I]) + '~' + pd.Series(st[J]), pd.Series(st[J]) + '~' + pd.Series(st[I]))
    df = pd.DataFrame(dict(i=I, j=J, sp=sp)); df = df[df.sp.isin(set(df.sp[ce[df.i.values] == ce[df.j.values]]))].reset_index(drop=True)
    print(kind, len(df), 'pairs', df.sp.nunique(), 'study pairs', flush=True)
    out = []
    for k, (i, j) in enumerate(zip(df.i.values, df.j.values)):
        b = fin[:, i] & fin[:, j]
        if b.sum() < 500: out.append((np.nan, np.nan)); continue
        xi, xj = X[b, i], X[b, j]; res = []
        for sg in [prim, -prim]:
            d = sg * xi > thr_d; r = sg * xj > thr_r
            if d.sum() < 10: res.append(np.nan); continue
            a = (d & r).sum() + .5; bb = (d & ~r).sum() + .5; c = (~d & r).sum() + .5; dd = (~d & ~r).sum() + .5
            res.append(np.log(a * dd / (bb * c)))
        out.append(tuple(res))
    df[['lor_primary', 'lor_reverse']] = np.array(out)
    units = E.unit.unique(); ucell = E.groupby('unit').cell.first().reindex(units).values; uidx = pd.Index(units).get_indexer(E.unit)
    def stat(cell, col):
        sc = cell[df.i.values] == cell[df.j.values]
        g = pd.DataFrame(dict(sp=df.sp, sc=sc, v=df[col])).dropna().groupby(['sp', 'sc']).v.mean().unstack()
        if True not in g: return np.nan, None
        g = g.dropna(); return (g[True] - g[False]).mean(), g
    rng = np.random.default_rng(2); res = []
    for col in ['lor_primary', 'lor_reverse']:
        obs, g = stat(ce, col)
        null = [stat(rng.permutation(ucell)[uidx], col)[0] for _ in range(nperm)]; null = np.array([x for x in null if not np.isnan(x)])
        stl = np.unique(np.concatenate([x.split('~') for x in g.index])); a = np.array([x.split('~')[0] for x in g.index]); b = np.array([x.split('~')[1] for x in g.index])
        dif = (g[True] - g[False]).values; bs = []
        for _ in range(1000):
            w = pd.Series(rng.choice(stl, len(stl))).value_counts(); ww = w.reindex(a).fillna(0).values * w.reindex(b).fillna(0).values
            if ww.sum() > 0: bs.append(np.average(dif, weights=ww))
        res.append(dict(kind=kind, dir=col.split('_')[1], lor_same=g[True].mean(), lor_diff=g[False].mean(), delta=obs, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                        perm_p=(np.sum(null >= obs) + 1) / (len(null) + 1), null_mean=null.mean(), n_studypairs=len(g), studypairs=';'.join(g.index)))
    df.to_csv(CA + f'sitelevel_pairs_{kind}.tsv.gz', sep='\t', index=False)
    return res
R2 = pd.DataFrame(site_level('inh') + site_level('stim'))
R2.to_csv(OUT + 'T3_H1_sitelevel_samecell.tsv', sep='\t', index=False); print(R2.drop(columns='studypairs').round(3).to_string())
