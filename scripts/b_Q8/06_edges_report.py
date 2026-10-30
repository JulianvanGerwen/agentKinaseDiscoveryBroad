"""b_Q8 step 6: edge-level preservation and rewiring meta-analysis (main 6 cohorts; OV/UCEC reported separately)."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from common import *
from scipy.stats import norm, spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
pd.set_option('display.width', 250)
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = np.empty_like(p); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
E = pd.read_csv(TQ + 'T3_compendium_edges.tsv', sep='\t'); E['rel'] = E.a + '|' + E.b; E = E.set_index('rel')
res = {}
glob_rows = []
for cond in ['padj_pc3', 'padj', 'padj_comp_pc3']:
    D = pd.read_parquet(CQ + f'edge_stats_{cond}.parquet')
    D['label'] = np.where(D.kind == 'site', 'fp_site', D.rel.map(E.label)); D['sign_c'] = np.where(D.kind == 'site', 1.0, np.sign(D.rel.map(E.r_meta)))
    D['r_comp'] = D.rel.map(E.r_meta)
    # ---- global preservation per cohort
    for coh, d in D.groupby('cohort'):
        for lab in ['edge', 'nonedge', 'fp_site']:
            x = d[d.label == lab]
            g = dict(cond=cond, cohort=coh, label=lab, n_rel=len(x),
                     signconc_T=(np.sign(x.r_T) == x.sign_c).mean() if lab != 'nonedge' else np.nan,
                     signconc_N=(np.sign(x.r_N) == x.sign_c).mean() if lab != 'nonedge' else np.nan,
                     pres_T=(x.sign_c * x.r_T).mean() if lab != 'nonedge' else x.r_T.abs().mean(),
                     pres_N=(x.sign_c * x.r_N).mean() if lab != 'nonedge' else x.r_N.abs().mean())
            if lab == 'edge':
                g['rho_comp_T'] = spearmanr(x.r_comp, x.r_T, nan_policy='omit')[0]; g['rho_comp_N'] = spearmanr(x.r_comp, x.r_N, nan_policy='omit')[0]
                g['rho_T_N'] = spearmanr(x.r_T, x.r_N, nan_policy='omit')[0]
            if lab == 'fp_site': g['rho_T_N'] = spearmanr(x.r_T, x.r_N, nan_policy='omit')[0]
            glob_rows.append(g)
    GP = pd.read_parquet(CQ + f'edge_global_{cond}.parquet')
    for coh, g in GP.groupby('cohort'):
        o = g[g.what == 'obs'].iloc[0]; pm = g[g.what == 'perm']
        for k in ['dpres_pairs', 'dpres_sites']:
            glob_rows.append(dict(cond=cond, cohort=coh, label='perm_' + k, pres_T=o[k], pres_N=np.nan,
                                  z_perm=(o[k] - pm[k].mean()) / pm[k].std(), p_perm=(1 + (np.abs(pm[k] - pm[k].mean()) >= abs(o[k] - pm[k].mean())).sum()) / (1 + len(pm))))
    # ---- per-relationship meta (main cohorts)
    M = D[D.cohort.isin(MAINC) & D.label.isin(['edge', 'fp_site'])].copy()
    M['sz'] = M.sign_c * M.z_perm; M['sdr'] = M.sign_c * (M.r_T - M.r_N)
    g = M.groupby('rel')
    R = pd.DataFrame(dict(kind=g.kind.first(), label=g.label.first(), k=g.sz.count(), r_comp=g.r_comp.first(),
                          rT_med=g.apply(lambda x: (x.sign_c * x.r_T).median()), rN_med=g.apply(lambda x: (x.sign_c * x.r_N).median()),
                          dr_med=g.sdr.median(), meta_z=g.sz.sum() / np.sqrt(g.sz.count()),
                          meta_z_fisher=g.apply(lambda x: (x.sign_c * x.z_fisher).sum() / np.sqrt(len(x))),
                          n_up=g.apply(lambda x: ((x.p_perm < 0.05) & (x.sz > 0)).sum()), n_down=g.apply(lambda x: ((x.p_perm < 0.05) & (x.sz < 0)).sum())))
    R = R[R.k >= 4]; R['meta_p'] = 2 * norm.sf(R.meta_z.abs()); R['meta_q'] = np.nan
    for lab in ['edge', 'fp_site']:
        m = R.label == lab; R.loc[m, 'meta_q'] = bh(R.loc[m, 'meta_p'])
    res[cond] = R
G = pd.DataFrame(glob_rows); G.to_csv(TQ + 'T7_edge_global_preservation.tsv', sep='\t', index=False)
print(G[G.label.isin(['edge', 'nonedge', 'fp_site'])].groupby(['cond', 'label'])[['signconc_T', 'signconc_N', 'pres_T', 'pres_N', 'rho_comp_T', 'rho_comp_N', 'rho_T_N']].median().round(3))
print(G[G.label.str.startswith('perm')].pivot_table(index=['cohort'], columns=['cond', 'label'], values='z_perm').round(2))
print(G[(G.cond == 'padj_pc3') & (G.label == 'edge')][['cohort', 'n_rel', 'signconc_T', 'signconc_N', 'pres_T', 'pres_N', 'rho_comp_T', 'rho_comp_N']].round(3).to_string())
R = res['padj_pc3'].copy()
R['meta_z_padj'] = res['padj'].meta_z; R['meta_z_comp_pc3'] = res['padj_comp_pc3'].meta_z
R['rewired'] = (R.meta_q < 0.05) & ((R.n_up >= 3) | (R.n_down >= 3)) & (np.sign(R.meta_z_comp_pc3) == np.sign(R.meta_z)) & (R.meta_z_comp_pc3.abs() > 1.96)
R['direction'] = np.where(R.meta_z > 0, 'stronger in tumour', 'weaker in tumour')
R.sort_values('meta_z').to_csv(TQ + 'T8_edge_rewiring.tsv', sep='\t')
for lab in ['edge', 'fp_site']:
    x = R[R.label == lab]
    print(lab, 'n', len(x), 'q<.05', (x.meta_q < .05).sum(), 'rewired(robust,>=3 coh)', x.rewired.sum(), 'up', (x.rewired & (x.meta_z > 0)).sum(), 'down', (x.rewired & (x.meta_z < 0)).sum())
    print(x[x.rewired].sort_values('meta_z')[['r_comp', 'rN_med', 'rT_med', 'dr_med', 'meta_z', 'n_up', 'n_down', 'meta_q', 'meta_z_padj', 'meta_z_comp_pc3']].round(3).head(40).to_string())
    print(x[x.rewired].sort_values('meta_z')[['r_comp', 'rN_med', 'rT_med', 'dr_med', 'meta_z', 'n_up', 'n_down', 'meta_q', 'meta_z_comp_pc3']].round(3).tail(25).to_string())
# rewired site edges by footprint
x = R[(R.label == 'fp_site')]; x = x.assign(fp=[r.split('>')[0] for r in x.index])
t = x.groupby('fp').agg(n=('meta_z', 'size'), n_rew_up=('rewired', lambda v: (v & (x.loc[v.index, 'meta_z'] > 0)).sum()),
                        n_rew_down=('rewired', lambda v: (v & (x.loc[v.index, 'meta_z'] < 0)).sum()), mean_meta_z=('meta_z', 'mean'),
                        rN=('rN_med', 'median'), rT=('rT_med', 'median'))
t.sort_values('mean_meta_z').to_csv(TQ + 'T9_site_rewiring_by_footprint.tsv', sep='\t'); print(t.sort_values('mean_meta_z').round(3).to_string())
