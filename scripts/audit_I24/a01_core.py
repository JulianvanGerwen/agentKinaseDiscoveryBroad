"""1) Independent recompute of I24 gains + IV decomposition; 2) instrument strength/composition, amplification exclusion,
per-receptor; correlated-error (bgPY) and abundance confounding; immune/stromal/SFK adjustment."""
from common_a import *
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
ST = {}; USED = []
for v in ['tc', 'tc_padj']:
    for c in COH:
        D, used = build(c, v)
        # amplification / receptor protein
        gi = cl.matrix(c, 'cnv_gistic'); pr = cl.matrix(c, 'protein')
        for g in ['EGFR', 'ERBB2', 'MET']:
            D[g + '_amp'] = (gi.loc[g].reindex(D.case_id).values == 2).astype(float) if g in gi.index else 0.
            D[g + '_prot'] = pr.loc[ACC[g]].reindex(D.index).values if ACC[g] in pr.index else np.nan
        for g in ['SHC1', 'PTPN11', 'GAB1', 'MAPK1', 'MAPK3']:
            D[g + '_prot'] = pr.loc[ACC[g]].reindex(D.index).values if ACC[g] in pr.index else np.nan
        D['ADAPT_prot'] = D[['SHC1_prot', 'PTPN11_prot', 'GAB1_prot']].mean(1); D['ERK_prot'] = D[['MAPK1_prot', 'MAPK3_prot']].mean(1)
        D['RTK_prot'] = D[['EGFR_prot', 'ERBB2_prot', 'MET_prot']].mean(1)
        D['PTPRC_prot'] = pr.loc['P08575'].reindex(D.index).values if 'P08575' in pr.index else np.nan
        D['cohort'] = c; D['variant'] = v
        ST[(v, c)] = D
        if v == 'tc':
            for k, s in used.items(): USED.append(dict(cohort=c, stage=k, sites=';'.join(s), n_sites=len(s), cover=D[k].notna().mean()))
A = pd.concat(ST.values()); A.to_parquet(f'{OUT}/tables/stages_audit.parquet')
U_ = pd.DataFrame(USED); U_.to_csv(f'{OUT}/tables/A0_stage_sites.tsv', sep='\t', index=False)
print(U_[U_.stage.isin(['RTK', 'SFK', 'bgPY'])].to_string())

def run(name, v, sub=lambda D: D, **kw):
    rows = []
    for c in COH:
        D = sub(ST[(v, c)].copy())
        r = decomp(D, **kw)
        if r is None: continue
        r.update(cohort=c, analysis=name, variant=v); rows.append(r)
    return rows

rows = []
for v in ['tc', 'tc_padj']:
    rows += run('main', v)
    rows += run('main_ERKSUB_na', v, y='ERK', cov=('purity', 'prolif'))  # placeholder identical (sanity)
    rows += run('no_amp', v, sub=lambda D: D[(D.EGFR_amp + D.ERBB2_amp + D.MET_amp) == 0])
    rows += run('EGFR_only', v, x='EGFRpY')
    rows += run('EGFR_only_no_amp', v, x='EGFRpY', sub=lambda D: D[D.EGFR_amp == 0])
    rows += run('adj_bgPY', v, cov=('purity', 'prolif', 'bgPY'))
    rows += run('adj_immune_stroma', v, cov=('purity', 'prolif', 'Immune', 'Stromal'))
    rows += run('adj_SFK', v, cov=('purity', 'prolif', 'SFK'))
    rows += run('adj_all', v, cov=('purity', 'prolif', 'bgPY', 'Immune', 'Stromal', 'SFK'))
    rows += run('outcome_ERK2', v, y='ERK2'); rows += run('outcome_ERK3', v, y='ERK3')
    rows += run('mediator_SHC1', v, m='ADAPT_SHC1'); rows += run('mediator_PTPN11_GAB1', v, m='ADAPT_PG')
rows += run('tc_adj_protein', 'tc', cov=('purity', 'prolif', 'ADAPT_prot', 'ERK_prot', 'RTK_prot'))
rows += run('tc_instr_RTKprotein', 'tc', x='RTK_prot', sub=lambda D: D[D.RTK.notna()])
R = pd.DataFrame(rows); R = R[R.analysis != 'main_ERKSUB_na']
R.to_csv(f'{OUT}/tables/A1_decomp_by_cohort.tsv', sep='\t', index=False)
M = meta_table(R, ['analysis', 'variant']); M.to_csv(f'{OUT}/tables/A1_decomp_meta.tsv', sep='\t', index=False)
cols = ['analysis', 'variant', 'k', 'n', 'medF', 'minF', 'fs', 'rf', 'iv', 'iv_lo', 'iv_hi', 'other', 'other_lo', 'other_hi', 'diff', 'diff_p', 'total']
print(M[cols].round(3).to_string())
print(R[R.analysis == 'main'][['variant', 'cohort', 'n', 'F', 'fs', 'rf', 'iv', 'iv_se', 'other', 'other_se', 'diff', 'diff_p_boot', 'total']].round(3).to_string())

# simple gains (I24 T2 analogue, OLS standardized partial slopes) for transitivity product check
g = []
for v in ['tc', 'tc_padj']:
    for c in COH:
        D = ST[(v, c)]
        for x, y in [('RTK', 'ADAPT'), ('ADAPT', 'ERK'), ('RTK', 'ERK')]:
            d = D[[x, y, 'purity', 'prolif']].dropna()
            if len(d) < 30: continue
            Z = d.apply(zs); X_ = np.c_[np.ones(len(d)), Z[[x, 'purity', 'prolif']].values]
            b = np.linalg.lstsq(X_, Z[y].values, rcond=None)[0]; r = Z[y].values - X_ @ b
            se = np.sqrt(r @ r / (len(d) - 4) * np.linalg.inv(X_.T @ X_)[1, 1])
            g.append(dict(variant=v, cohort=c, edge=f'{x}>{y}', b=b[1], se=se, n=len(d)))
G = pd.DataFrame(g); G.to_csv(f'{OUT}/tables/A1_gains_by_cohort.tsv', sep='\t', index=False)
gm = [dict(variant=v, edge=e, **meta_re(d.b, d.se)) for (v, e), d in G.groupby(['variant', 'edge'])]
gm = pd.DataFrame(gm); gm.to_csv(f'{OUT}/tables/A1_gains_meta.tsv', sep='\t', index=False); print(gm.round(3).to_string())
# same 7 cohorts for ADAPT>ERK
gm7 = [dict(variant=v, edge=e, **meta_re(d.b, d.se)) for (v, e), d in G[G.cohort.isin(R[R.analysis == 'main'].cohort.unique())].groupby(['variant', 'edge'])]
print('7-cohort gains'); print(pd.DataFrame(gm7).round(3).to_string())
