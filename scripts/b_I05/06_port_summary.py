"""Summaries/contrasts of portability pairs (05). Study-pair weighting + study-cluster bootstrap (as audit_I19 boot_diff)."""
import numpy as np, pandas as pd
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
R = pd.read_csv(CA + 'portability_pairs.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
rng = np.random.default_rng(3)


def sp_mean(d, col='lor'):
    return d.groupby(['si', 'sj'])[col].mean()


def boot_contrast(d1, d2=None, nb=1000):
    """mean over study pairs of (d1 - d2) (or d1) ; cluster bootstrap over studies"""
    a = sp_mean(d1)
    if d2 is not None:
        b = sp_mean(d2); x = (a - b).dropna()
    else: x = a.dropna()
    if len(x) == 0: return np.nan, np.nan, np.nan, 0
    si = x.index.get_level_values(0).values; sj = x.index.get_level_values(1).values
    st = np.unique(np.r_[si, sj]); bs = []
    for _ in range(nb):
        w = pd.Series(rng.choice(st, len(st))).value_counts()
        ww = w.reindex(si).fillna(0).values * w.reindex(sj).fillna(0).values
        if ww.sum() > 0: bs.append(np.average(x.values, weights=ww))
    return x.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(x)


out = []
for tdef, tims in [('timing', ['early_sus', 'early_trans', 'late', 'EARLY']), ('timing1', ['early', 'late'])]:
    for kind in ['inh', 'stim', 'all']:
        d0 = R[(R.tdef == tdef) & ((R.kind == kind) | (kind == 'all'))].copy()
        if tdef == 'timing':
            e = d0[d0.timing.str.startswith('early')].copy(); e['timing'] = 'EARLY'
            # pool early_sus and early_trans: n_disc-weighted LOR within pair
            e = e.groupby(['series', 'si', 'sj', 'sj', 'prim', 'timing', 'bin', 'dirv'], as_index=False).apply(
                lambda g: pd.Series(dict(lor=np.average(g.lor, weights=g.n_disc), rep=np.average(g.rep, weights=g.n_disc), base=g.base.iloc[0], n_disc=g.n_disc.sum())), include_groups=False) if len(e) else e
            d0 = pd.concat([d0, e], ignore_index=True)
        for b in ['0.5-1', '1-2', '2-99', '>1']:
            db = d0[d0.bin.isin(['1-2', '2-99'])] if b == '>1' else d0[d0.bin == b]
            for pr in ['primary', 'reverse']:
                for tm in tims:
                    g = db[(db.prim == pr) & (db.timing == tm)]
                    m, lo, hi, n = boot_contrast(g)
                    out.append(dict(tdef=tdef, kind=kind, bin=b, contrast=f'{pr}|{tm}', est=m, lo=lo, hi=hi, n_studypairs=n,
                                    rep=sp_mean(g, 'rep').mean() if len(g) else np.nan, base=sp_mean(g, 'base').mean() if len(g) else np.nan))
            E_ = 'EARLY' if tdef == 'timing' else 'early'
            C = {'primary: early-late': (('primary', E_), ('primary', 'late')), 'reverse: early-late': (('reverse', E_), ('reverse', 'late')),
                 'early: primary-reverse': (('primary', E_), ('reverse', E_)), 'late: primary-reverse': (('primary', 'late'), ('reverse', 'late')),
                 'primary-early minus reverse-late': (('primary', E_), ('reverse', 'late'))}
            if tdef == 'timing':
                C['primary: transient-late'] = (('primary', 'early_trans'), ('primary', 'late'))
                C['reverse: transient-late'] = (('reverse', 'early_trans'), ('reverse', 'late'))
                C['primary: sustained-transient'] = (('primary', 'early_sus'), ('primary', 'early_trans'))
            for name, (c1, c2) in C.items():
                g1 = db[(db.prim == c1[0]) & (db.timing == c1[1])]; g2 = db[(db.prim == c2[0]) & (db.timing == c2[1])]
                # within-series contrast: match on series x replication exp
                k = ['series', 'sj', 'bin', 'dirv'] if False else ['series', 'si', 'sj']
                m, lo, hi, n = boot_contrast(g1, g2)
                out.append(dict(tdef=tdef, kind=kind, bin=b, contrast='D ' + name, est=m, lo=lo, hi=hi, n_studypairs=n))
O = pd.DataFrame(out)
O.to_csv(OUT + 'T3_portability_timing.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
print(O[(O.tdef == 'timing') & O.bin.isin(['1-2', '2-99', '>1'])].round(3).to_string())
print(O[(O.tdef == 'timing1') & O.bin.isin(['>1'])].round(3).to_string())
