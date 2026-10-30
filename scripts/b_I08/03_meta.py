#!/usr/bin/env python3
"""I08 step 3: meta-analysis of coupling (Q1), null controls (Q2), classification, reporter sites (Q4).
Inputs: cache/coupling_long.parquet, joint_r2.parquet. Outputs: tables/*.tsv"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *
from scipy import stats

L = pd.read_parquet(f'{OUT}/cache/coupling_long.parquet')
kin = pd.read_csv(f'{R}/model_corev3/kinases.tsv', sep='\t'); fam = dict(zip(kin.kinase, kin.manning_family)); grp = dict(zip(kin.kinase, kin.manning_group))
canon = pd.read_csv(f'{OUT}/tables/canonical_site.tsv', sep='\t').set_index('kinase').canonical_site.to_dict()
ann = pd.read_csv(f'{OUT}/tables/kinase_site_annotation.tsv', sep='\t').set_index('site_key')
L['site'] = L.feature.str.split('|').str[0].where(L.feature.str.contains(r'\|'))
L['mode'] = L.feature.str.split('|').str[1]
L['ftype'] = np.where(L.feature.isin(['prot', 'rna']), L.feature,
                      np.where(L.site == L.kinase.map(canon), 'canon_' + L['mode'].fillna(''), 'site_' + L['mode'].fillna('')))


def meta_block(g):
    out = {}
    for col, kk in [('r', 0), ('r_partial', 3)]:
        m = meta(g[col].values, g.n.values, kk)
        sfx = '' if col == 'r' else '_p'
        out.update({f'r{sfx}': m['r'], f'lo{sfx}': m['lo'], f'hi{sfx}': m['hi'], f'p{sfx}': m['p'], f'I2{sfx}': m['I2']})
        # leave-one-cohort-out: number of LOCO fits that stay p<0.05 with the same sign
        if m['k'] >= 3:
            ok = 0
            for c in g.cohort:
                gg = g[g.cohort != c]; mm = meta(gg[col].values, gg.n.values, kk)
                ok += (mm['p'] < 0.05) and np.sign(mm['r']) == np.sign(m['r'])
            out[f'loco_ok{sfx}'] = ok
        rr = g[col].values; n = g.n.values
        z = np.arctanh(np.clip(rr, -.999, .999)) * np.sqrt(np.maximum(n - 3 - kk, 1))
        out[f'n_rep{sfx}'] = int(((np.abs(z) > 1.96) & (np.sign(rr) == np.sign(m['r']))).sum())
    out['k'] = len(g)
    if 'null_sd' in g and g.null_sd.notna().any():
        zp = (g.r - g.null_mean) / g.null_sd; zpp = (g.r_partial - g.nullp_mean) / g.nullp_sd
        out['perm_z'] = zp.sum() / np.sqrt(zp.notna().sum()); out['perm_z_p'] = zpp.sum() / np.sqrt(zpp.notna().sum())
    if 'pct_other' in g and g.pct_other.notna().any():
        out['pct_other'] = g.pct_other.median(); out['pct_other_p'] = g.pctp_other.median()
        out['best_other_mode'] = g.best_otherp.mode().iat[0]; out['other_median_r'] = g.other_median_r.median()
    return pd.Series(out)


res = []
for (src, k, f), g in L.groupby(['source', 'kinase', 'feature']):
    if g.n.ge(20).sum() < 3: continue
    s = meta_block(g); s['source'] = src; s['kinase'] = k; s['feature'] = f; s['ftype'] = g.ftype.iat[0]
    res.append(s)
M = pd.DataFrame(res)
M['fdr'] = np.nan; M['fdr_p'] = np.nan
for (src, ft), idx in M.groupby(['source', 'ftype']).groups.items():
    M.loc[idx, 'fdr'] = bh(M.loc[idx, 'p']); M.loc[idx, 'fdr_p'] = bh(M.loc[idx, 'p_p'])
M['family'] = M.kinase.map(fam); M['group'] = M.kinase.map(grp)
M.to_csv(f'{OUT}/tables/meta_all.tsv', sep='\t', index=False)
print(M.groupby(['source', 'ftype']).size())
