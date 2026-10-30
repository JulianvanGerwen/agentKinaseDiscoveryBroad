"""H4: variance in cross-study reverse-layer discordance (per pair LOR_up, |log2FC|>1) explained by |delta axis shift| between the two experiments."""
import numpy as np, pandas as pd
OUT = '/workspace/results/b_I25/tables/'
R = pd.read_csv('/workspace/cache/b_I25/pairs_all.tsv', sep='\t', dtype={'si': str, 'sj': str})
E = pd.read_csv('/workspace/cache/b_I25/E_test.tsv', sep='\t', index_col=0)
F = pd.read_csv(OUT + 'T3_H3_exp_features.tsv', sep='\t', index_col='exp').reindex(E.exp.values).reset_index(drop=True)
R = R[R.bin == '1-99']
P = R.pivot_table(index=['i', 'j', 'si', 'sj', 'group'], columns='dir', values='lor').reset_index().dropna(subset=['up'])
# symmetrise: average (i,j),(j,i)
P['a'] = P[['i', 'j']].min(1); P['b'] = P[['i', 'j']].max(1)
P = P.groupby(['a', 'b']).agg(up=('up', 'mean'), down=('down', 'mean'), si=('si', 'first'), sj=('sj', 'first'), group=('group', 'first')).reset_index()
FE = ['shift_M', 'shift_MKI', 'shift_INTER', 'shift_MET24', 'shift_MET30', 'shift_CULT3D', 'shift_SERUM', 'shift_OIS', 'shift_HIPPO', 'shift_AMPKmotif', 'frac_down', 'median', 'P_prim', 'R_rev']
for f in FE: P['d_' + f] = (F[f].values[P.a] - F[f].values[P.b]).__abs__()
P['same_cell'] = (E.cell.values[P.a] == E.cell.values[P.b]).astype(float)
P['d_time'] = np.abs(np.log10(E.time_min.fillna(60).values[P.a] + 1) - np.log10(E.time_min.fillna(60).values[P.b] + 1))
P.to_csv(OUT + 'T4_H4_pairs.tsv', sep='\t', index=False)
studies = np.unique(np.r_[P.si, P.sj]); rng = np.random.default_rng(11)
def r2(d, cols, y='up', base=()):
    d = d.dropna(subset=list(cols) + [y] + list(base))
    def fit(cs):
        X = np.c_[np.ones(len(d)), d[list(cs)].values] if cs else np.ones((len(d), 1))
        b = np.linalg.lstsq(X, d[y].values, rcond=None)[0]; res = d[y].values - X @ b
        return 1 - res.var() / d[y].var()
    return fit(list(base) + list(cols)) - (fit(list(base)) if base else 0), len(d)
def boot(cols, base=(), y='up', nb=500):
    out = []
    for _ in range(nb):
        w = pd.Series(rng.choice(studies, len(studies))).value_counts()
        ww = (w.reindex(P.si).fillna(0).values * w.reindex(P.sj).fillna(0).values).astype(int)
        d = P.loc[np.repeat(P.index.values, ww)]
        if len(d) > 20: out.append(r2(d, cols, y, base)[0])
    return np.percentile(out, [2.5, 97.5])
rows = []
for f in FE + ['same_cell', 'd_time']:
    c = f if f in ('same_cell', 'd_time') else 'd_' + f
    for y, base in [('up', ()), ('up', ('down',)), ('down', ())]:
        v, n = r2(P, [c], y, base); lo, hi = boot([c], base, y, 300)
        d = P.dropna(subset=[c, y]); rho = np.corrcoef(d[c], d[y])[0, 1]
        rows.append(dict(predictor=c, outcome=y, adjusted_for='+'.join(base), R2=v, lo=lo, hi=hi, r=rho, n_pairs=n))
cc = ['d_' + f for f in FE[:10]] + ['d_frac_down', 'd_median']
for y, base in [('up', ()), ('up', ('down',))]:
    v, n = r2(P, cc, y, base); lo, hi = boot(cc, base, y, 300)
    rows.append(dict(predictor='JOINT_all_axes', outcome=y, adjusted_for='+'.join(base), R2=v, lo=lo, hi=hi, n_pairs=n))
v, n = r2(P, ['down'], 'up'); rows.append(dict(predictor='down_LOR(same pair)', outcome='up', R2=v, lo=boot(['down'], (), 'up', 300)[0], hi=np.nan, n_pairs=n))
T = pd.DataFrame(rows); T.to_csv(OUT + 'T4_H4_discordance_R2.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(len(P), 'pairs', P.si.nunique(), studies.size, 'studies'); print(T.round(3).to_string())
