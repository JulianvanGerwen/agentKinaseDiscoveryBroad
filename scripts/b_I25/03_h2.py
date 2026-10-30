"""H2 (+H4 site-level): cross-study same-class replication LOR for up (reverse) vs down (primary) sites, after removing axis sets
vs matched random removal (same residue x coverage x variance strata composition)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I25'); from common_i25 import *
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import pair_rows
OUT = '/workspace/results/b_I25/tables/'
A, E = load_test('inh'); AX = axes(A.index); Av = A.values.astype(np.float32)
BINS = [(1, 2), (1, 99)]
cls_groups = {c: (E.cls == c).values for c in E.cls.unique() if E[E.cls == c].study.nunique() >= 2}
cell_groups = {f'{c}|{k}': ((E.cls == c) & (E.cell == k)).values for (c, k), d in E.groupby(['cls', 'cell']) if d.study.nunique() >= 2}
strat = AX.stratum.fillna('NA').values
rng = np.random.default_rng(7)
def matched_random(setmask, base):
    out = np.zeros(len(setmask), bool)
    s = pd.Series(np.arange(len(setmask)))
    cand = base & ~setmask
    for st, n in pd.Series(strat[setmask & base]).value_counts().items():
        pool = np.where(cand & (strat == st))[0]
        if len(pool): out[rng.choice(pool, min(n, len(pool)), replace=False)] = True
    return out
def run(mask, groups):
    R = pair_rows(Av, E, groups, bins=BINS, site_mask=mask, null_pool=None, min_both=200, min_disc=5)
    return R
def summ(R):
    sp = R.groupby(['bin', 'dir', 'si', 'sj']).lor.mean().groupby(['bin', 'dir']).mean().unstack('dir')
    sp['gap'] = sp['down'] - sp['up']; return sp
def boot(R, nb=300):
    st = np.unique(np.r_[R.si, R.sj]); spm = R.groupby(['bin', 'dir', 'si', 'sj']).lor.mean().reset_index(); res = []
    for _ in range(nb):
        w = pd.Series(rng.choice(st, len(st))).value_counts()
        ww = w.reindex(spm.si).fillna(0).values * w.reindex(spm.sj).fillna(0).values
        d = spm.assign(w=ww)[ww > 0]
        res.append(d.groupby(['bin', 'dir']).apply(lambda g: np.average(g.lor, weights=g.w), include_groups=False))
    return pd.concat(res, axis=1)
SETS = ['CC_union', 'CC_Mup', 'CC_Mdown', 'CC_MKI', 'CC_inter', 'S_MET24', 'S_MET30', 'S_AMPKmotif', 'S_CULT3D', 'S_HIPPO', 'S_SERUM', 'S_OIS']
rows = []; allR = {}
for gname, groups, base in [('class_xstudy_UCC', cls_groups, AX.U_CC.values), ('class_xstudy_all', cls_groups, np.ones(len(AX), bool)),
                            ('samecell_xstudy_UCC', cell_groups, AX.U_CC.values)]:
    R0 = run(base, groups); s0 = summ(R0); allR[(gname, 'none')] = R0
    b0 = boot(R0)
    for bn in s0.index:
        rows.append(dict(scope=gname, removed='none', bin=bn, up=s0.loc[bn, 'up'], down=s0.loc[bn, 'down'], gap=s0.loc[bn, 'gap'],
                         up_lo=np.percentile(b0.loc[(bn, 'up')], 2.5), up_hi=np.percentile(b0.loc[(bn, 'up')], 97.5),
                         gap_lo=np.percentile(b0.loc[(bn, 'down')] - b0.loc[(bn, 'up')], 2.5), gap_hi=np.percentile(b0.loc[(bn, 'down')] - b0.loc[(bn, 'up')], 97.5), n_removed=0))
    sets = SETS if gname == 'class_xstudy_all' else ['CC_union', 'CC_Mup', 'CC_MKI']
    for sname in sets:
        sm = AX[sname].values
        R1 = run(base & ~sm, groups); s1 = summ(R1); allR[(gname, sname)] = R1
        Rin = run(base & sm, groups); sin = summ(Rin) if len(Rin) else None
        rnd = []
        for k in range(20 if gname != 'samecell_xstudy_UCC' else 10):
            rm = matched_random(sm, base); rnd.append(summ(run(base & ~rm, groups)))
        for bn in s1.index:
            ru = np.array([r.loc[bn, 'up'] for r in rnd]); rg = np.array([r.loc[bn, 'gap'] for r in rnd])
            rows.append(dict(scope=gname, removed=sname, bin=bn, up=s1.loc[bn, 'up'], down=s1.loc[bn, 'down'], gap=s1.loc[bn, 'gap'], n_removed=int((base & sm).sum()),
                             delta_up=s1.loc[bn, 'up'] - s0.loc[bn, 'up'], gap_closure=1 - s1.loc[bn, 'gap'] / s0.loc[bn, 'gap'],
                             rnd_delta_up_mean=ru.mean() - s0.loc[bn, 'up'], rnd_delta_up_hi=np.percentile(ru, 97.5) - s0.loc[bn, 'up'],
                             rnd_gap_closure_mean=1 - rg.mean() / s0.loc[bn, 'gap'],
                             within_up=sin.loc[bn, 'up'] if sin is not None and bn in sin.index else np.nan,
                             within_down=sin.loc[bn, 'down'] if sin is not None and bn in sin.index else np.nan))
        print(gname, sname, 'done', flush=True)
T = pd.DataFrame(rows); T.to_csv(OUT + 'T2_H2_removal.tsv', sep='\t', index=False)
allR[('class_xstudy_all', 'none')].to_csv('/workspace/cache/b_I25/pairs_all.tsv', sep='\t', index=False)
E.to_csv('/workspace/cache/b_I25/E_test.tsv', sep='\t')
pd.set_option('display.width', 250); print(T.round(3).to_string())
