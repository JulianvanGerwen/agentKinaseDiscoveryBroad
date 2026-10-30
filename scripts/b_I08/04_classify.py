#!/usr/bin/env python3
"""I08 step 4: per-kinase classification (Q1/Q2). Uses meta_all.tsv (partial = adjusted for purity, ESTIMATE, prolif).
abundance  : own protein partial r>0, FDR<0.05, LOCO-stable in >=80% of LOCO fits, beats random-substrate null (perm z>3, meanz)
             and beats other kinases' protein (median percentile >= 0.75)
PTM site   : any PSP-activating / A-loop site (protein-residualised) partial r>0, FDR<0.05, LOCO >=80%
classes    : abundance / abundance+PTM / PTM / uncoupled (+ 'no site measured')
variability: SNR = SD(observed activity) / median SD(200 random same-size substrate sets), per cohort, median (meanz).
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *
M = pd.read_csv(f'{OUT}/tables/meta_all.tsv', sep='\t')
ann = pd.read_csv(f'{OUT}/tables/kinase_site_annotation.tsv', sep='\t').set_index('site_key')
M['site'] = M.feature.str.split('|').str[0].where(M.feature.str.contains(r'\|'))
M['act_site'] = M.site.map(lambda s: bool(s in ann.index and (ann.loc[s, 'psp_activating'] or ann.loc[s, 'in_aloop'])) if isinstance(s, str) else False)

# SNR from nulls (meanz_tc_padj)
act = pd.read_parquet(f'{OUT}/cache/act_meanz_tc_padj.parquet'); snr = []
for c in cl.COHORTS:
    nz = np.load(f'{OUT}/cache/null_tc_padj_{c}.npz', allow_pickle=True); cols = list(nz['samples'])
    for i, k in enumerate(nz['kin']):
        a = act.loc[k, cols].values.astype(float); nsd = np.nanmedian(np.nanstd(nz['null'][i], 1))
        snr.append((k, c, np.nanstd(a) / nsd))
snr = pd.DataFrame(snr, columns=['kinase', 'cohort', 'snr']).groupby('kinase').snr.median()

J = pd.read_parquet(f'{OUT}/cache/joint_r2.parquet')
out = []
for src in ['meanz_tc_padj', 'corev3_tc_padj', 'meanz_tc', 'corev3_tc']:
    S = M[M.source == src]
    P = S[S.ftype == 'prot'].set_index('kinase'); Rn = S[S.ftype == 'rna'].set_index('kinase')
    A = S[S.act_site & S.feature.str.endswith('|padj')].copy()
    A['fdr_act'] = bh(A.p_p)
    A = A.sort_values('r_p', ascending=False).drop_duplicates('kinase').set_index('kinase')
    ks = sorted(set(P.index) | set(A.index) | set(Rn.index))
    for k in ks:
        d = dict(source=src, kinase=k, family=S[S.kinase == k].family.iat[0], group=S[S.kinase == k].group.iat[0])
        for nm, T in [('prot', P), ('rna', Rn)]:
            if k in T.index:
                t = T.loc[k]
                d.update({f'{nm}_r': t.r, f'{nm}_rp': t.r_p, f'{nm}_lo': t.lo_p, f'{nm}_hi': t.hi_p, f'{nm}_fdr': t.fdr_p,
                          f'{nm}_loco': t.loco_ok_p / t.k, f'{nm}_nrep': t.n_rep_p, f'{nm}_k': t.k, f'{nm}_permz': t.get('perm_z_p', np.nan),
                          f'{nm}_permz_raw': t.get('perm_z', np.nan)})
                if nm == 'prot': d.update(prot_pct_other=t.pct_other_p, prot_pct_other_raw=t.pct_other, prot_best_other=t.best_other_mode)
        if k in A.index:
            t = A.loc[k]
            d.update(site=t.site, site_r=t.r, site_rp=t.r_p, site_lo=t.lo_p, site_hi=t.hi_p, site_fdr=t.fdr_act,
                     site_loco=t.loco_ok_p / t.k, site_nrep=t.n_rep_p, site_k=t.k, site_permz=t.get('perm_z_p', np.nan))
        j = J[(J.source == src) & (J.kinase == k)]
        d['R2adj_rna_prot_site_mean'] = j.R2_adj.mean() if len(j) else np.nan
        out.append(d)
C = pd.DataFrame(out)
meanz = C.source.str.startswith('meanz')
ab = (C.prot_fdr < .05) & (C.prot_rp > 0) & (C.prot_loco >= .8) & (C.prot_pct_other >= .75) & (~meanz | (C.prot_permz > 3))
pt = (C.site_fdr < .05) & (C.site_rp > 0) & (C.site_loco >= .8) & (~meanz | (C.site_permz > 3))
C['abund'] = ab; C['ptm'] = pt
C['cls'] = np.select([ab & pt, ab, pt, C.site.isna()], ['abundance+PTM', 'abundance', 'PTM', 'uncoupled(no site)'], 'uncoupled')
C['snr'] = C.kinase.map(snr)
C.to_csv(f'{OUT}/tables/kinase_classification_all_sources.tsv', sep='\t', index=False)
# consensus over the two methods on the primary input (tc_padj)
W = C[C.source.isin(['meanz_tc_padj', 'corev3_tc_padj'])].pivot(index='kinase', columns='source', values='cls')
W.columns = ['cls_corev3', 'cls_meanz']
W['consensus'] = np.where(W.cls_meanz.isna(), 'corev3 only: ' + W.cls_corev3.astype(str),
                          np.where(W.cls_meanz == W.cls_corev3, W.cls_meanz, 'discordant'))
W.to_csv(f'{OUT}/tables/kinase_class_consensus.tsv', sep='\t')
for s in ['meanz_tc_padj', 'corev3_tc_padj', 'meanz_tc', 'corev3_tc']:
    print(s, C[C.source == s].cls.value_counts().to_dict())
print(W[W.cls_meanz.notna()].consensus.value_counts())
pd.set_option('display.width', 250)
print(W[W.cls_meanz.notna() & (W.cls_meanz != 'uncoupled(no site)')].sort_values('consensus').to_string())
