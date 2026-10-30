"""Extra contrasts for task 2 (sustained primary vs reverse; same-cell excluded) and summary of task 3a (motif-matched vs not).
Fast study-cluster bootstrap (numpy)."""
import numpy as np, pandas as pd
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
rng = np.random.default_rng(4)


def boot(x, nb=1000):
    """x: Series indexed by (si, sj) -> mean, CI via study cluster bootstrap"""
    x = x.dropna()
    if len(x) < 3: return x.mean() if len(x) else np.nan, np.nan, np.nan, len(x)
    si = x.index.get_level_values(0).values; sj = x.index.get_level_values(1).values
    st, inv = np.unique(np.r_[si, sj], return_inverse=True); ii, jj = inv[:len(x)], inv[len(x):]
    bs = []
    for _ in range(nb):
        cnt = np.bincount(rng.integers(0, len(st), len(st)), minlength=len(st)); w = cnt[ii] * cnt[jj]
        if w.sum() > 0: bs.append(np.average(x.values, weights=w))
    return x.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(x)


sp = lambda d: d.groupby(['si', 'sj']).lor.mean()
# ---- task 2 extra
R = pd.read_csv(CA + 'portability_pairs.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
R = R[(R.tdef == 'timing') & ~R.same_cell]
out = []
for kind in ['inh', 'stim', 'all']:
    d0 = R if kind == 'all' else R[R.kind == kind]
    for b in [['1-2'], ['2-99'], ['1-2', '2-99']]:
        d = d0[d0.bin.isin(b)]
        g = lambda p, t: d[(d.prim == p) & (d.timing == t)]
        for t in ['early_sus', 'early_trans', 'late']:
            for p in ['primary', 'reverse']:
                out.append(dict(kind=kind, bin='+'.join(b), contrast=f'{p}|{t}', **dict(zip(['est', 'lo', 'hi', 'n'], boot(sp(g(p, t)))))))
            out.append(dict(kind=kind, bin='+'.join(b), contrast=f'D {t}: primary-reverse', **dict(zip(['est', 'lo', 'hi', 'n'], boot(sp(g('primary', t)) - sp(g('reverse', t)))))))
        out.append(dict(kind=kind, bin='+'.join(b), contrast='D primary: sustained-transient', **dict(zip(['est', 'lo', 'hi', 'n'], boot(sp(g('primary', 'early_sus')) - sp(g('primary', 'early_trans')))))))
        out.append(dict(kind=kind, bin='+'.join(b), contrast='D reverse: sustained-transient', **dict(zip(['est', 'lo', 'hi', 'n'], boot(sp(g('reverse', 'early_sus')) - sp(g('reverse', 'early_trans')))))))
        out.append(dict(kind=kind, bin='+'.join(b), contrast='D primary: transient-late', **dict(zip(['est', 'lo', 'hi', 'n'], boot(sp(g('primary', 'early_trans')) - sp(g('primary', 'late')))))))
        out.append(dict(kind=kind, bin='+'.join(b), contrast='DD (prim-rev)[sus] - (prim-rev)[late]', **dict(zip(['est', 'lo', 'hi', 'n'], boot(
            sp(g('primary', 'early_sus')) - sp(g('reverse', 'early_sus')) - sp(g('primary', 'late')) + sp(g('reverse', 'late')))))))
O = pd.DataFrame(out); O.to_csv(OUT + 'T3b_portability_persistence_diffcell.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
print(O[O.bin == '1-2+2-99'].round(2).to_string())
# ---- task 3a
D = pd.read_csv(CA + 'direct_portability_pairs.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
D = D[~D.same_cell]
out = []
for kind in ['inh', 'stim']:
    for split in ['motif', 'psp']:
        for p in ['primary', 'reverse']:
            for b in [['0.5-1'], ['1-2'], ['2-99'], ['1-2', '2-99']]:
                d = D[(D.kind == kind) & (D.split == split) & (D.prim == p) & D.bin.isin(b)]
                m1 = sp(d[d.matched]); m0 = sp(d[~d.matched])
                out.append(dict(kind=kind, split=split, prim=p, bin='+'.join(b), lor_matched=m1.mean(), lor_unmatched=m0.mean(),
                                **dict(zip(['D_est', 'lo', 'hi', 'n'], boot(m1 - m0)))))
        # interaction: (matched-unmatched)[primary] - (matched-unmatched)[reverse]
        d = D[(D.kind == kind) & (D.split == split) & D.bin.isin(['1-2', '2-99'])]
        x = (sp(d[d.matched & (d.prim == 'primary')]) - sp(d[~d.matched & (d.prim == 'primary')])) - (sp(d[d.matched & (d.prim == 'reverse')]) - sp(d[~d.matched & (d.prim == 'reverse')]))
        out.append(dict(kind=kind, split=split, prim='interaction prim-rev', bin='1-2+2-99', **dict(zip(['D_est', 'lo', 'hi', 'n'], boot(x)))))
O3 = pd.DataFrame(out); O3.to_csv(OUT + 'T4_direct_portability.tsv', sep='\t', index=False)
print(O3.round(2).to_string())
# per class (motif, primary, >1)
pc = []
for (kind, c), d in D[(D.split == 'motif') & D.bin.isin(['1-2', '2-99'])].groupby(['kind', 'cls']):
    for p in ['primary', 'reverse']:
        dd = d[d.prim == p]; m1 = sp(dd[dd.matched]); m0 = sp(dd[~dd.matched])
        x = (m1 - m0).dropna()
        pc.append(dict(kind=kind, cls=c, prim=p, D=x.mean(), n_sp=len(x), n_studies=len(set(x.index.get_level_values(0)) | set(x.index.get_level_values(1))),
                       n_matched=int(dd[dd.matched].n_disc.sum())))
PC = pd.DataFrame(pc).pivot_table(index=['kind', 'cls'], columns='prim', values=['D', 'n_sp']).round(2)
PC.to_csv(OUT + 'T4b_direct_portability_by_class.tsv', sep='\t'); print(PC.to_string())
# ---- placebo kinase sets (07b): same D (matched - unmatched) per draw
Pz = pd.read_csv(CA + 'direct_portability_PLACEBO_pairs.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
Pz = Pz[~Pz.same_cell]
pl = []
for (kind, p, dr), d in Pz.groupby(['kind', 'prim', 'draw']):
    x = (sp(d[d.matched]) - sp(d[~d.matched])).dropna()
    pl.append(dict(kind=kind, prim=p, draw=dr, D=x.mean(), n=len(x)))
PL = pd.DataFrame(pl); PL.to_csv(OUT + 'T4c_placebo_kinase_sets.tsv', sep='\t', index=False)
print(PL.groupby(['kind', 'prim']).D.agg(['mean', 'min', 'max']).round(2))
