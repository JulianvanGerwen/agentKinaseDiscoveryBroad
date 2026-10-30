#!/usr/bin/env python3
"""I08 step 5 (Q4): activity-reporter phosphosites on kinases.
For every phosphosite on kinase k (protein-residualised 'padj' and raw), meta partial r with k's substrate-inferred
activity (from meta_all.tsv) + specificity: per cohort, correlation of the site with the activity of EVERY kinase ->
percentile of the own kinase (and best non-family kinase). Reporter = FDR<0.05 (BH over all kinase sites), LOCO>=80%,
>=4 cohorts, positive, specificity percentile >= 0.9, in BOTH meanz and corev3 (tc_padj).
Output: tables/reporter_sites.tsv, tables/site_specificity.tsv"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *
kin = pd.read_csv(f'{R}/model_corev3/kinases.tsv', sep='\t'); k2acc = dict(zip(kin.kinase, kin.uniprot)); fam = dict(zip(kin.kinase, kin.manning_family))
acc2k = dict(zip(kin.uniprot, kin.kinase))
ann = pd.read_csv(f'{OUT}/tables/kinase_site_annotation.tsv', sep='\t').set_index('site_key')
M = pd.read_csv(f'{OUT}/tables/meta_all.tsv', sep='\t')
M['site'] = M.feature.str.split('|').str[0].where(M.feature.str.contains(r'\|')); M['mode'] = M.feature.str.split('|').str[1]
canon = pd.read_csv(f'{OUT}/tables/canonical_site.tsv', sep='\t').set_index('kinase').canonical_site.to_dict()

spec = []
for src in ['meanz_tc_padj', 'corev3_tc_padj']:
    meth, var = src.split('_', 1)
    ACT = pd.read_parquet(f'{OUT}/cache/act_{meth}_{var}.parquet')
    for c in cl.COHORTS:
        T = [t for t in cl.tumour_ids(c) if t in ACT.columns]
        C = covariates(c, T); Cr = np.column_stack([rank(C[x]) for x in C])
        A = ACT[T]; A = A[A.notna().sum(1) >= 20]
        Ar = np.vstack([resid(rank(r), Cr) for r in A.values])
        kn = np.array(A.index)
        for mode in ['padj']:
            X = cl.matrix(c, 'input_tc_padj').reindex(columns=T)
            X = X[X.index.str.split('_').str[0].isin(acc2k)]
            for s, x in zip(X.index, X.values.astype(float)):
                k = acc2k[s.split('_')[0]]
                if k not in kn or np.isfinite(x).sum() < 20: continue
                xp = resid(rank(x), Cr); m = np.isfinite(xp)
                Am = Ar[:, m]; ok = np.isfinite(Am).all(1)
                rr = np.full(len(kn), np.nan)
                xx = xp[m] - xp[m].mean(); AA = Am[ok] - Am[ok].mean(1, keepdims=True)
                rr[ok] = AA @ xx / (np.sqrt((AA ** 2).sum(1)) * np.sqrt((xx ** 2).sum()) + 1e-12)
                for i in np.where(~ok)[0]:          # activity rows with missing samples: pairwise-complete
                    rr[i] = pcor(Ar[i], xp)[0]
                own = kn == k; other = ~own & np.isfinite(rr)
                nonfam = other & (np.array([fam.get(q) for q in kn]) != fam.get(k))
                j = np.nanargmax(np.where(nonfam, rr, -9))
                spec.append(dict(source=src, cohort=c, site=s, kinase=k, mode=mode, r_own=rr[own][0],
                                 pct=np.mean(rr[other] < rr[own][0]), pct_nonfam=np.mean(rr[nonfam] < rr[own][0]),
                                 best_nonfam=kn[j], best_nonfam_r=rr[j]))
        print(src, c, len(spec), flush=True)
S = pd.DataFrame(spec)
S.to_csv(f'{OUT}/tables/site_specificity_long.tsv', sep='\t', index=False)
Sg = S.groupby(['source', 'site', 'mode']).agg(spec_pct=('pct', 'median'), spec_pct_nonfam=('pct_nonfam', 'median'),
                                               best_nonfam=('best_nonfam', lambda x: x.mode().iat[0]), best_nonfam_r=('best_nonfam_r', 'median')).reset_index()
X = M[M.site.notna() & M.source.isin(['meanz_tc_padj', 'corev3_tc_padj'])].copy()
X = X.merge(Sg, on=['source', 'site', 'mode'], how='left')
X['fdr_all'] = np.nan
for (src, md), idx in X.groupby(['source', 'mode']).groups.items(): X.loc[idx, 'fdr_all'] = bh(X.loc[idx, 'p_p'])
X = X.join(ann[['psp_activating', 'psp_inhibitory', 'psp_regulatory', 'in_aloop', 'psp_function', 'res', 'pos']], on='site')
X['canonical'] = X.site == X.kinase.map(canon)
X['loco_frac'] = X.loco_ok_p / X.k
X['pass'] = (X.fdr_all < .05) & (X.r_p > 0) & (X.loco_frac >= .8) & (X.k >= 4) & (X.spec_pct >= .9) & \
            (~X.source.str.startswith('meanz') | (X.perm_z_p > 3))
X.to_csv(f'{OUT}/tables/kinase_site_activity_coupling.tsv', sep='\t', index=False)
P = X[X['mode'] == 'padj']
W = P.pivot_table(index=['kinase', 'site'], columns='source', values=['r_p', 'lo_p', 'hi_p', 'fdr_all', 'k', 'n_rep_p', 'spec_pct', 'pass', 'perm_z_p'], aggfunc='first')
W.columns = [f'{a}__{b.split("_")[0]}' for a, b in W.columns]
W = W.reset_index().join(ann[['psp_activating', 'psp_inhibitory', 'psp_regulatory', 'in_aloop', 'psp_function']], on='site')
W['canonical_site'] = W.kinase.map(canon)
cr = P[P.canonical].set_index(['source', 'kinase']).r_p
W['canon_rp__meanz'] = [cr.get(('meanz_tc_padj', k), np.nan) for k in W.kinase]
W['canon_rp__corev3'] = [cr.get(('corev3_tc_padj', k), np.nan) for k in W.kinase]
W['both_pass'] = (W['pass__meanz'] == True) & (W['pass__corev3'] == True)
W['one_pass_other_na'] = ((W['pass__meanz'] == True) & W['pass__corev3'].isna()) | ((W['pass__corev3'] == True) & W['pass__meanz'].isna())
W['status'] = np.where(W.psp_activating | W.in_aloop, 'known activating/A-loop', np.where(W.psp_regulatory, 'PSP regulatory (other)', 'novel (no PSP regulatory annotation)'))
W.sort_values('r_p__corev3', ascending=False).to_csv(f'{OUT}/tables/reporter_sites_wide.tsv', sep='\t', index=False)
R_ = W[W.both_pass | W.one_pass_other_na].copy()
R_['beats_canonical'] = [(not (r.site == r.canonical_site)) and (np.nanmax([r['r_p__meanz'], r['r_p__corev3']]) > np.nanmax([r.canon_rp__meanz, r.canon_rp__corev3, -1])) for _, r in R_.iterrows()]
R_.to_csv(f'{OUT}/tables/reporter_sites.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 30)
print(R_[['kinase', 'site', 'r_p__meanz', 'r_p__corev3', 'k__corev3', 'spec_pct__meanz', 'spec_pct__corev3', 'both_pass', 'canonical_site', 'canon_rp__corev3', 'beats_canonical', 'status']].round(2).to_string())
