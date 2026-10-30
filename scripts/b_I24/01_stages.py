"""Build per-tumour stage scores (mean of cohort-centred log2 sites) for tc and tc_padj, plus covariates."""
from common import *
psp = pd.read_csv('/workspace/data/Kinase_Substrate_Dataset', sep='\t', skiprows=3, encoding='latin-1')
psp = psp[(psp.KIN_ORGANISM == 'human') & (psp.SUB_ORGANISM == 'human') & psp.KIN_ACC_ID.isin(['P28482', 'P27361'])]
excl_genes = set(U) | {'SOS1', 'SOS2', 'RAF1', 'BRAF', 'ARAF', 'MAP2K1', 'MAP2K2', 'DUSP4', 'DUSP6', 'DUSP5', 'DUSP1', 'MKNK1', 'MKNK2',
                       'RPS6KA2', 'RPS6KA4', 'RPS6KA5', 'RPS6KA6', 'MAPK1', 'MAPK3', 'SPRY2', 'SPRY4', 'SPRED1', 'SPRED2', 'GAB2', 'SHC1'}
excl_acc = set(U.values())
psp = psp[~psp.SUB_GENE.isin(excl_genes) & ~psp.SUB_ACC_ID.isin(excl_acc) & ~psp.SUB_ACC_ID.str.contains('-')]
erk_sites = sorted(set(psp.SUB_ACC_ID + '_' + psp.SUB_MOD_RSD))
print('PSP ERK substrate sites (excl cascade proteins):', len(erk_sites))
stage_sites = set(s for v in STAGES.values() for s in v)
assert not (set(erk_sites) & stage_sites)
cov_rows, cover = [], []
for variant in ['tc', 'tc_padj']:
    allrows = []
    for c in COHORTS:
        T = cl.tumour_ids(c)
        X = cl.matrix(c, f'input_{variant}')
        X = X[[t for t in T if t in X.columns]]
        df = pd.DataFrame(index=X.columns)
        for st, sites in STAGES.items():
            ss = [s for s in sites if s in X.index]
            if not ss:
                df[st] = np.nan; df[st + '_n'] = 0; continue
            df[st] = X.loc[ss].mean(0, skipna=True); df[st + '_n'] = X.loc[ss].notna().sum(0)
            if variant == 'tc_padj': cover.append((c, st, len(ss), ';'.join(s for s in ss), df[st].notna().mean()))
        es = [s for s in erk_sites if s in X.index]
        Z = X.loc[es]; Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1), axis=0)
        n = Z.notna().sum(0)
        df['ERKSUB'] = Z.mean(0).where(n >= 5); df['ERKSUB_n'] = n
        if variant == 'tc_padj': cover.append((c, 'ERKSUB', len(es), 'PSP ERK1/2 substrates', df['ERKSUB'].notna().mean()))
        df['cohort'] = c; df['case_id'] = [s[:-2] for s in df.index]
        allrows.append(df)
    A = pd.concat(allrows)
    # covariates
    cl_ = cl.clinical().drop_duplicates('case_id').set_index('case_id')
    ann = pd.concat([cl.table(c, 'annotations') for c in COHORTS]).drop_duplicates('case_id').set_index('case_id')
    A['purity'] = cl_.reindex(A.case_id)['WES_purity'].values
    A['prolif'] = ann.reindex(A.case_id)['HALLMARK_E2F_TARGETS'].values
    A['estimate'] = ann.reindex(A.case_id)['ESTIMATE_ESTIMATEScore'].values
    for c in COHORTS:  # fill missing purity with cohort median
        m = A.cohort == c
        for v in ['purity', 'prolif', 'estimate']:
            A.loc[m, v] = A.loc[m, v].fillna(A.loc[m, v].median())
    A.to_parquet(f'{OUT}/tables/stages_{variant}.parquet')
    print(variant, A.shape)
cv = pd.DataFrame(cover, columns=['cohort', 'stage', 'n_sites_measured', 'sites', 'frac_tumours'])
cv.to_csv(f'{OUT}/tables/T1_stage_coverage.tsv', sep='\t', index=False)
print(cv.pivot(index='stage', columns='cohort', values='frac_tumours').round(2).to_string())
