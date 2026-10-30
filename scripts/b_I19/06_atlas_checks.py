"""Checks on atlas hits: (1) substrate drivers (which sites carry the signal; curated vs predicted edges),
(2) normalisation artefact (global median / fraction-down of each experiment; uncentred raw substrate log2FC),
(3) RAF paradox by BRAF status, (4) curated-only re-score."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import *
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/b_I19/"
E = experiments()
X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]
Xc = center(X, E)
G = pd.read_csv(CA + 'exp_globalstats.tsv', sep='\t', index_col=0).reindex(E.exp)
n = X.notna().sum(); G['frac_down1'] = (X < -1).sum() / n; G['frac_up1'] = (X > 1).sum() / n; G['median'] = X.median()
G.to_csv(CA + 'exp_globalstats.tsv', sep='\t')
ed = pd.read_csv('/workspace/scripts/infra_compendium/kinact_mvp_pkg/model/edges.tsv.gz', sep='\t')
R = pd.read_csv(OUT + 'tables/T2_compensation_atlas_all.tsv', sep='\t')
c = R[R.method == 'corev3'].merge(R[R.method == 'psp_z'][['cls', 'kinase', 'est', 'z', 'p_null_match_up']], on=['cls', 'kinase'], how='left', suffixes=('', '_psp'))
hits = c[(c.q < 0.05) & (c.est > 0) & (c.k_studies >= 2)].copy()
# add positive-control pairs for comparison
ctrl = c[c.set_index(['cls', 'kinase']).index.isin([('AKT', 'MAPK1'), ('AKT', 'MAPK3'), ('rapalog', 'AKT1'), ('MEK', 'AKT1')])]
hits = pd.concat([hits, ctrl])
A = pd.read_parquet(CA + 'act_corev3_exp.parquet').reindex(columns=E.exp)
Ac = center(A, E)
out = []
for r in hits.itertuples():
    eK = E[(E.cls == r.cls) & ~E.bset.apply(lambda b: r.kinase in b)].exp
    subs = ed[ed.kinase == r.kinase]
    s_in = subs[subs.site.isin(X.index)]
    v = Xc.loc[s_in.site, eK]
    mv = v.mean(1); nmeas = v.notna().sum(1)
    top = mv[nmeas >= max(2, len(eK) // 3)].sort_values(ascending=False)
    drivers = ';'.join(f'{s}({ed_src}):{m:.2f}' for s, m, ed_src in [(s, top[s], s_in.set_index('site').source.get(s)) for s in top.index[:6]])
    cur = s_in[s_in.source == 'cur'].site
    vc = Xc.loc[cur, eK]; cur_mean = np.nanmean(vc.values) if vc.notna().sum().sum() else np.nan
    frac_pos_sub = (top > 0).mean()
    # artefact: corr of kinase activity with global median shift / frac down across ALL experiments; and within class
    a = Ac.loc[r.kinase] if r.kinase in Ac.index else pd.Series(np.nan, index=E.exp)
    ok = a.notna()
    rho_down = pd.Series(a[ok].values).corr(pd.Series(G.frac_down1[ok].values), method='spearman')
    rho_med = pd.Series(a[ok].values).corr(pd.Series(G['median'][ok].values), method='spearman')
    # raw (uncentred, un-z-scored) substrate log2FC under K minus global median of the same experiment
    raw = X.loc[s_in.site, eK]; raw_mean = np.nanmean(raw.values); medK = G.loc[eK, 'median'].mean(); fdK = G.loc[eK, 'frac_down1'].mean()
    out.append(dict(cls=r.cls, kinase=r.kinase, est=r.est, z=r.z, k=r.k_studies, n=r.n_exp, est_psp=r.est_psp, n_sub_edges=len(subs), n_sub_meas=int((nmeas > 0).sum()),
                    n_curated_meas=int(vc.notna().any(axis=1).sum()), curated_sub_mean_centered=cur_mean, frac_substrates_up=frac_pos_sub,
                    raw_sub_lfc_mean=raw_mean, mean_exp_median_shift=medK, mean_frac_down=fdK, rho_activity_vs_fracdown=rho_down, rho_activity_vs_median=rho_med, drivers=drivers))
H = pd.DataFrame(out)
H.to_csv(OUT + 'tables/T2_hits_checks.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_colwidth', 160)
print(H.round(2).to_string())
# global artefact test: number of up kinases (z>2) per experiment vs fraction down
Az = A.where(pd.read_parquet(CA + 'act_corev3_exp.parquet').notna())
nup = (Ac > 2).sum(); ndn = (Ac < -2).sum()
from scipy import stats
print('n_up_kinases vs frac_down1 spearman', stats.spearmanr(nup, G.frac_down1), ' vs median', stats.spearmanr(nup, G['median']))
print('n_up vs n_down kinases', stats.spearmanr(nup, ndn))
# RAF paradox by BRAF status (judgement: V600 lines)
v600 = {'COLO205', 'WIDR', 'A375', 'WM239A', 'HT29', 'WM266'}
raf = E[E.cls == 'RAF'].copy(); raf['braf'] = np.where(raf.cell_line.isin(v600), 'V600', 'WT')
for k in ['MAPK1', 'MAPK3', 'RPS6KA1', 'RPS6KA3']:
    a = Ac.loc[k]; raf[k] = a[raf.exp].values
P = pd.read_parquet(CA + 'act_psp_z_exp.parquet').reindex(columns=E.exp); Pc = center(P, E)
for k in ['MAPK1', 'MAPK3']:
    raf[k + '_psp'] = Pc.loc[k][raf.exp].values
raf[['exp', 'braf', 'time_min', 'MAPK1', 'MAPK3', 'RPS6KA1', 'RPS6KA3', 'MAPK1_psp', 'MAPK3_psp']].to_csv(OUT + 'tables/T2_raf_paradox.tsv', sep='\t', index=False)
print(raf[['exp', 'braf', 'time_min', 'MAPK1', 'MAPK3', 'RPS6KA1', 'MAPK1_psp', 'MAPK3_psp']].round(2).to_string())
print(raf.groupby('braf')[['MAPK1', 'MAPK3', 'RPS6KA1', 'MAPK1_psp']].median())
