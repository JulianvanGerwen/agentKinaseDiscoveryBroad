"""Empirical null for the motif filter on the discovery set: swap each kinase's motif for a random kinase of a different
Manning group but the same residue class (ST vs Y); compare counts of motif>=90 hits under nested robustness filters."""
import numpy as np, pandas as pd
from lib import *
rng = np.random.default_rng(2)
M = pd.read_parquet(f'{RES}/tables/T3_meta_kinase_site.parquet')
motif = pd.read_parquet(f'{CACHE}/motif_pct.parquet')
kl = pd.read_csv('<HACKATHON_ROOT>/kinaseSubMain/results/labels/kinase_list.tsv', sep='\t')
grp = kl.drop_duplicates('gene').set_index('gene').manning_group.to_dict()
isYk = motif.loc[[i for i in motif.index if '_Y' in i]].notna().any() & ~motif.loc[[i for i in motif.index if '_Y' not in i][:5000]].notna().any()
base = (M.k >= 2) & (M.q_meta < 0.05) & (M.beta_iv > 0) & (M.n_pos >= 2) & ~M.cis & ~M.psp_edge & M.motif_pct.notna()
M['nb_robust'] = (M.n_nb_av == 0) | (M.n_nb_ok >= np.ceil(M.n_nb_av / 2))
filters = {'base': base, '+nb_robust': base & M.nb_robust, '+nb_robust+focal_t>1': base & M.nb_robust & (M.mean_t_fo > 1),
           '+nb+focal+naive_t>1': base & M.nb_robust & (M.mean_t_fo > 1) & (M.mean_t_nv > 1)}
rows = []
for nm, f in filters.items():
    sub = M[f]
    obs = (sub.motif_pct >= 90).sum()
    nulls = []
    for r in range(50):
        sw = {}
        for k in sub.kgene.unique():
            c = [x for x in motif.columns if grp.get(x, 'a') != grp.get(k, 'b') and isYk.get(x) == isYk.get(k)]
            sw[k] = rng.choice(c)
        v = np.array([motif.at[s, sw[k]] for s, k in zip(sub.site, sub.kgene)])
        nulls.append(np.nansum(v >= 90) / np.isfinite(v).sum() * len(sub))
    rows.append(dict(filter=nm, n=len(sub), own_motif90=obs, own_frac=obs / max(len(sub), 1), null_mean=np.mean(nulls),
                     null_sd=np.std(nulls), excess=obs - np.mean(nulls), est_FDR_motif=min(1, np.mean(nulls) / max(obs, 1)),
                     exp_frac_uniform=0.10))
r = pd.DataFrame(rows); print(r.round(3).to_string())
r.to_csv(f'{RES}/tables/T4b_motif_null.tsv', sep='\t', index=False)
hc = M[filters['+nb_robust+focal_t>1'] & (M.motif_pct >= 90)].sort_values('p_meta')
hc.to_csv(f'{RES}/tables/T4c_high_confidence_edges.tsv', sep='\t', index=False)
print(hc[['kgene', 'site', 'sgene', 'k', 'n_pos', 'cohorts_pos', 'beta_iv', 'se_iv', 'z', 'q_meta', 'motif_pct', 'psp_anysite', 'mean_t_nv', 'mean_t_fo']].head(30).to_string())
