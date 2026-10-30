"""H2/H4: does baseline expression (DepMap RNA; CCLE proteome; HPA RNA fallback) or Chronos dependency of kinase k in cell c predict the
reverse-layer motif enrichment of k in units (study x cell) of that cell?
Unit profile = mean of experiment profiles (dpct; centred per experiment), then centred per kinase across units (within kind).
Global test: per unit Spearman over kinases of profile vs baseline z (z per gene across the cells used); mean over units; null = permute
cell->baseline assignment (2000). Within-study test: profiles and baseline z centred within study (studies with >=2 mapped cells).
Kinase-family version: baseline of motif column k = max baseline z among kinases whose Cantley percentile vectors correlate r>=0.8 with k.
Targeted: pre-declared reverse kinases, across cells (one value per cell), Spearman."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
CA = '/workspace/cache/b_I22/'; OUT = '/workspace/results/b_I22/tables/'
P = pd.read_parquet(CA + 'profiles.parquet'); M = pd.read_csv(CA + 'profiles_meta.tsv', sep='\t', dtype={'study': str})
kst = list(P.columns)
B = {'rna': pd.read_parquet(CA + 'rna_depmap.parquet'), 'prot': pd.read_parquet(CA + 'prot_ccle.parquet'),
     'hpa': pd.read_parquet(CA + 'rna_hpa.parquet'), 'chronos': pd.read_parquet(CA + 'chronos.parquet')}
MP = pd.read_parquet('/workspace/results/b_I05/cache/motif_pct.parquet', columns=kst).sample(30000, random_state=0)
Cm = np.corrcoef(MP.values.T.astype(float)); fam = {k: [kst[j] for j in np.where(Cm[i] >= 0.8)[0]] for i, k in enumerate(kst)}
pd.Series({k: ','.join(v) for k, v in fam.items()}).to_csv(OUT + 'motif_families_r0.8.tsv', sep='\t')
DECL = ['MAP4K4', 'TNIK', 'MINK1', 'MAP4K2', 'MAP4K3', 'MAP4K5', 'NEK1', 'NEK3', 'NEK5', 'NEK9', 'IRAK4', 'TTBK1', 'TTBK2', 'EIF2AK2', 'MAP3K2', 'MAP3K3', 'MARK1', 'MARK2', 'MARK3', 'MARK4', 'SIK1', 'SIK2', 'SIK3']

def unit_profiles(dr, defn='dpct', cen='raw'):
    ix = M.index[(M.defn == defn) & (M.cen == cen) & (M.dir == dr)].values
    m = M.loc[ix].reset_index(drop=True); Z = P.loc[ix].values.astype(float); Z -= Z.mean(1, keepdims=True)
    for k in m.kind.unique():
        s = (m.kind == k).values; Z[s] -= Z[s].mean(0)
    m['unit'] = m.kind + '|' + m.study + '|' + m.cell
    U = pd.DataFrame(Z, columns=kst).groupby(m.unit.values).mean()
    um = m.groupby('unit').agg(kind=('kind', 'first'), study=('study', 'first'), cell=('cell', 'first'), n=('exp', 'size')).loc[U.index]
    return U, um

def baseline_z(src, cells, family=False, within=None):
    Bm = B[src]; cells = [c for c in cells if c in Bm.index]
    g = [k for k in kst if k in Bm.columns]
    X = Bm.loc[cells, g].astype(float)
    X = X.loc[:, X.std() > 0]
    Zb = (X - X.mean()) / X.std()
    if family:
        Zb = pd.DataFrame({k: Zb[[x for x in fam[k] if x in Zb.columns]].max(1) for k in kst if any(x in Zb.columns for x in fam[k])})
    return Zb

def test(dr, src, family=False, within=False, sign=1, nperm=2000, defn='dpct', cen='raw'):
    U, um = unit_profiles(dr, defn, cen)
    um = um[um.cell.isin(B[src].index)]
    if within:
        ok = um.groupby('study').cell.transform('nunique') >= 2; um = um[ok]
    U = U.loc[um.index]
    cells = sorted(um.cell.unique()); Zb = baseline_z(src, cells, family)
    ks = [k for k in kst if k in Zb.columns]; Uv = U[ks].copy()
    if within:  # centre profiles within study across units
        Uv = Uv - Uv.groupby(um.study.values).transform('mean')
    Ur = np.apply_along_axis(rankdata, 1, Uv.values)
    def stat(assign):  # assign: dict cell->row of Zb
        rs = []
        for u, (st, c) in enumerate(zip(um.study, um.cell)):
            b = Zb.loc[assign[c], ks].values.astype(float)
            if within:
                cs = um.cell[um.study == st].unique(); b = b - np.mean([Zb.loc[assign[x], ks].values for x in cs], 0)
            f = np.isfinite(b)
            if f.sum() < 20: continue
            rs.append(np.corrcoef(rankdata(Uv.values[u][f]), rankdata(b[f]))[0, 1])
        return np.nanmean(rs) * sign, rs
    obs, rs = stat({c: c for c in cells})
    rng = np.random.default_rng(3); null = []
    for _ in range(nperm):
        perm = rng.permutation(cells)
        if within:  # permute only among cells within each study is too restrictive for 2-3 cells; permute globally
            pass
        null.append(stat(dict(zip(cells, perm)))[0])
    null = np.array(null)
    return dict(dir=dr, src=src, family=family, within=within, defn=defn, cen=cen, mean_rho=obs, null_mean=null.mean(), null_sd=null.std(),
                perm_p=(np.sum(null >= obs) + 1) / (nperm + 1), n_units=len(um), n_cells=len(cells), n_kinases=len(ks))

rows = []
for dr in ['reverse', 'primary']:
    for src in ['rna', 'prot', 'chronos']:
        for fam_ in [False, True]:
            for w in [False, True]:
                try:
                    r = test(dr, src, fam_, w, nperm=500)
                    rows.append(r); print(r, flush=True)
                except Exception as e:
                    print('fail', dr, src, fam_, w, e)
T = pd.DataFrame(rows); T.to_csv(OUT + 'T5_H2_H4_baseline_vs_profile.tsv', sep='\t', index=False)
pd.set_option('display.width', 220); print(T.round(3).to_string())
# targeted: one value per cell (mean over its units), Spearman across cells
U, um = unit_profiles('reverse')
cellU = U.groupby(um.cell.values).mean()
tr = []
for src in ['rna', 'prot', 'chronos']:
    Bm = B[src]
    for k in DECL:
        if k not in Bm.columns or k not in cellU.columns: continue
        cs = [c for c in cellU.index if c in Bm.index]
        if len(cs) < 6: continue
        rho, p = spearmanr(cellU.loc[cs, k], Bm.loc[cs, k], nan_policy='omit')
        famk = [x for x in fam[k] if x in Bm.columns]
        rho2, p2 = spearmanr(cellU.loc[cs, k], Bm.loc[cs, famk].max(1)) if famk else (np.nan, np.nan)
        tr.append(dict(src=src, kinase=k, n_cells=len(cs), rho=rho, p=p, rho_family_max=rho2, p_family=p2))
TR = pd.DataFrame(tr); TR.to_csv(OUT + 'T6_H2_targeted_kinases.tsv', sep='\t', index=False); print(TR.round(3).to_string())
