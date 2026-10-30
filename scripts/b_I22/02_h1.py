"""H1: is the reverse-layer motif profile a stable property of the cell line across unrelated perturbation classes?
Pairwise Pearson r of profiles (centred per experiment, then per kinase within kind x dir). Matched contrast:
within each study pair (s1<s2, different studies), mean r(same cell, diff class) - mean r(diff cell, diff class); averaged over study pairs that have both.
Null: block permutation of cell labels over study x cell units (1000). CI: study cluster bootstrap (1000).
Also within-study contrasts for decryptM (independent per-drug controls) and Hijazi (shared DMSO control; raw vs centred)."""
import numpy as np, pandas as pd, itertools, sys
CA = '/workspace/cache/b_I22/'; OUT = '/workspace/results/b_I22/tables/'
P = pd.read_parquet(CA + 'profiles.parquet'); M = pd.read_csv(CA + 'profiles_meta.tsv', sep='\t', dtype={'study': str})
Mask = pd.read_parquet(CA + 'measured_mask.parquet')

def prep(defn, cen, dr, kinds=('inh', 'stim'), kcentre=True):
    ix = M.index[(M.defn == defn) & (M.cen == cen) & (M.dir == dr) & M.kind.isin(kinds)].values
    m = M.loc[ix].reset_index(drop=True); Z = P.loc[ix].values.astype(float)
    Z = Z - Z.mean(1, keepdims=True)
    if kcentre:
        for k in m.kind.unique():
            s = (m.kind == k).values; Z[s] -= Z[s].mean(0)
    Z /= np.linalg.norm(Z, axis=1, keepdims=True)
    m['cls2'] = m.kind + ':' + m.cls; m['unit'] = m.study + '|' + m.cell
    return m, Z @ Z.T

def pairs(m, R):
    I, J = np.triu_indices(len(m), 1)
    d = pd.DataFrame(dict(i=I, j=J, r=R[I, J], ss=(m.study.values[I] == m.study.values[J]), dc=(m.cls2.values[I] != m.cls2.values[J]),
                          ki=m.kind.values[I], kj=m.kind.values[J]))
    s1 = np.where(m.study.values[I] < m.study.values[J], m.study.values[I], m.study.values[J]); s2 = np.where(m.study.values[I] < m.study.values[J], m.study.values[J], m.study.values[I])
    d['sp'] = pd.Series(s1) + '~' + pd.Series(s2)
    return d

def matched(d, cell, I, J, wstudy=False):
    sc = cell[I] == cell[J]
    g = pd.DataFrame(dict(sp=d.sp.values, sc=sc, r=d.r.values)).groupby(['sp', 'sc']).r.mean().unstack()
    if True not in g or False not in g: return np.nan, g.iloc[:0]
    g = g.dropna(); return (g[True] - g[False]).mean(), g

def run(defn, cen, dr, sel, label, kinds=('inh', 'stim'), nperm=1000, nboot=1000, within=False):
    m, R = prep(defn, cen, dr, kinds)
    d = pairs(m, R)
    d = d[sel(d)].reset_index(drop=True)
    if within: d['sp'] = m.study.values[d.i]  # within-study: group by study
    I, J = d.i.values, d.j.values
    cell = m.cell.values
    obs, g = matched(d, cell, I, J)
    if np.isnan(obs): return None
    # permutation: shuffle cell labels across units
    units = m.unit.unique(); ucell = m.groupby('unit').cell.first().reindex(units).values; uidx = pd.Index(units).get_indexer(m.unit)
    rng = np.random.default_rng(1); null = []
    for _ in range(nperm):
        pc = rng.permutation(ucell)[uidx]
        v, _ = matched(d, pc, I, J)
        if not np.isnan(v): null.append(v)
    null = np.array(null)
    # bootstrap over studies (study pairs weighted by product of study multiplicities)
    st = np.unique(np.concatenate([np.array(x.split('~')) for x in g.index]))
    diffs = (g[True] - g[False]); a = np.array([x.split('~')[0] for x in g.index]); b = np.array([x.split('~')[-1] for x in g.index])
    bs = []
    for _ in range(nboot):
        w = pd.Series(rng.choice(st, len(st))).value_counts()
        ww = w.reindex(a).fillna(0).values * (w.reindex(b).fillna(0).values if not within else 1)
        if ww.sum() > 0: bs.append(np.average(diffs.values, weights=ww))
    return dict(test=label, defn=defn, cen=cen, dir=dr, r_same=g[True].mean(), r_diff=g[False].mean(), delta=obs, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                perm_p=(np.sum(null >= obs) + 1) / (len(null) + 1), null_mean=null.mean(), n_units=len(units), n_studypairs=len(g),
                n_pairs_same=int((cell[I] == cell[J]).sum()))

rows = []
for defn in ['dpct', 'lor90']:
    for cen in ['raw', 'cen']:
        for dr in ['reverse', 'primary']:
            tests = [
                ('A cross-study, diff class (inh+stim)', lambda d: (~d.ss) & d.dc, ('inh', 'stim'), False),
                ('B cross-study, diff class, inh only', lambda d: (~d.ss) & d.dc, ('inh',), False),
                ('C cross-study, inh vs stim', lambda d: (~d.ss) & (d.ki != d.kj), ('inh', 'stim'), False),
                ('D cross-study, any class (inh+stim)', lambda d: (~d.ss), ('inh', 'stim'), False),
                ('E cross-study, diff class, stim only', lambda d: (~d.ss) & d.dc, ('stim',), False),
            ]
            for lab, sel, kinds, w in tests:
                r = run(defn, cen, dr, sel, lab, kinds)
                if r: rows.append(r)
            # within-study tests: decryptM (independent controls) and Hijazi (shared control)
            for st, lab in [('36926954', 'F within decryptM, diff class (indep. controls)'), ('31959955', 'G within Hijazi, diff class (shared DMSO)'),
                            ('31959955', 'H within Hijazi, SAME drug (shared DMSO)')]:
                if 'SAME' in lab:
                    m0, _ = prep(defn, cen, dr)
                    ag = m0.agent.values
                    sel = (lambda st_, ag_: (lambda d: d.ss & (d.sp.str.startswith(st_ + '~')) & (ag_[d.i.values] == ag_[d.j.values]) | d.ss & d.sp.str.startswith(st_ + '~') & d.dc))(st, ag)
                else:
                    sel = (lambda st_: (lambda d: d.ss & d.sp.str.startswith(st_ + '~') & d.dc))(st)
                r = run(defn, cen, dr, sel, lab, within=True, nboot=0 + 1)
                if r: r['lo'] = r['hi'] = np.nan; rows.append(r)
            print(defn, cen, dr, 'done', flush=True)
T = pd.DataFrame(rows); T.to_csv(OUT + 'T1_H1_cell_identity.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(T.round(3).to_string())
