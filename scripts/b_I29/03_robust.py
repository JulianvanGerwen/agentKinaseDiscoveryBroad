"""I29 task 1: robustness, peptide forms, proteins, within-EGFR DUSP6 analysis (LUAD)."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from lib import *
A = load('LUAD')
A['pTpY_minus_monoY'] = A.pERK_broad - A.pERK_monoY
A['MAPK1_pTpY_minus_pY'] = z(A.MAPK1_pTpY) - z(A.MAPK1_pY); A['MAPK3_pTpY_minus_pY'] = z(A.MAPK3_pTpY) - z(A.MAPK3_pY)
A['TP53'] = A.TP53_mut.astype(float)
B = ['purity', 'prolif']; rows = []
# ---- H1/H2/H3 + outcome panel, base model
OUTS = ['stage_RTK', 'stage_EGFRpY', 'stage_ADAPT', 'stage_ERK', 'pERK_broad', 'MAPK1_pTpY', 'MAPK3_pTpY', 'pERK_monoY', 'MAPK1_pY', 'MAPK3_pY',
        'pTpY_minus_monoY', 'MAPK1_pTpY_minus_pY', 'MAPK3_pTpY_minus_pY', 'MAPK1_pTpY_raw', 'MAPK3_pTpY_raw', 'stage_ERK_tc',
        'stage_RSK', 'stage_ERKSUB', 'stage_ERKSUB_tc', 'MAPKout_fp', 'fp_EGF', 'fp_MAPKaxis', 'stage_AKTSUB', 'stage_MTORC1',
        'prot_MAPK1', 'prot_MAPK3', 'MAPK1_protBroad', 'MAPK3_protBroad', 'prot_DUSP6', 'prot_DUSP4', 'prot_SPRY2', 'prot_SPRY4', 'prot_SPRED1', 'prot_SPRED2', 'prot_ETV4', 'prot_ETV5', 'prot_EGFR', 'DUSP6_protBroad', 'DUSP4_protBroad',
        'ETVc'] + [c for c in A.columns if c.startswith('rna_')] + ['PROGENy_MAPK', 'PROGENy_EGFR']
OUTS = [o for o in OUTS if o in A.columns and A[o].notna().sum() >= 20]
for o in OUTS:
    rows += row(ols(A, o, B), test='panel', outcome=o, model='base (vs WTall)')
    rows += row(ols(A, o, B, groups=('EGFR',), ref='RAS'), test='panel', outcome=o, model='EGFR vs RAS-path')
# ---- H4 covariate robustness on primaries + Broad pTpY + ERKSUB + ETVc
PRIM = ['stage_ERK', 'pERK_broad', 'stage_ERKSUB', 'ETVc']
A['EGFR_clean'] = A.EGFR_clean.fillna(False).astype(bool); A['WT_clean'] = A.WT_clean.fillna(False).astype(bool)
A['grp_clean'] = np.where(A.EGFR_clean, 'EGFR', np.where(A.WT_clean, 'WT', np.where((A.group == 'RAS'), 'RAS', 'drop')))
A['grp_allele'] = np.where(A.EGFR_allele == 'ex19del', 'ex19', np.where(A.EGFR_allele == 'L858R', 'L858R', np.where(A.group == 'EGFR', 'EGFRother', A.group)))
for o in PRIM:
    for m, cv in [('+smoking', B + ['smoker']), ('+stage', B + ['stage_num']), ('+sex', B + ['male']), ('+TP53', B + ['TP53']),
                  ('+all clinical (smoking,stage,sex,TP53)', B + ['smoker', 'stage_num', 'male', 'TP53']),
                  ('+all clinical +ISCH_noERK +log excision', B + ['smoker', 'stage_num', 'male', 'TP53', 'ISCH_noERK', 'log_excision_f'])]:
        rows += row(ols(A, o, cv), test='H4 covariates', outcome=o, model=m)
    rows += row(ols(A, o, B, gcol='grp_clean'), test='H4 clean groups', outcome=o, model='EGFR ex19/L858R no RAS vs WT no fusion/EGFRamp/ERBB2')
    rows += row(ols(A, o, B + ['smoker'], gcol='grp_clean'), test='H4 clean groups', outcome=o, model='clean + smoking')
    rows += row(ols(A, o, B, groups=('ex19', 'L858R', 'EGFRother', 'RAS'), gcol='grp_allele'), test='H4 allele', outcome=o, model='by allele vs WTall')
    rows += row(ols(A, o, B, groups=('ex19',), ref='L858R', gcol='grp_allele'), test='H4 allele', outcome=o, model='ex19 vs L858R')
    rows += row(ols(A, o, B, mask=A.SMOKING_STATUS == 'Non-Smoker'), test='H4 strata', outcome=o, model='never-smokers only')
    rows += row(ols(A, o, B, mask=A.SMOKING_STATUS == 'Smoker'), test='H4 strata', outcome=o, model='smokers only')
    # KRAS only (exclude NF1/BRAF/NRAS) vs WT
    A['grp_kras'] = np.where(A.KRAS_mut == 1, 'KRAS', A.group)
    rows += row(ols(A, o, B, groups=('KRAS',), gcol='grp_kras', mask=A.EGFR_mut == 0), test='H3 KRAS', outcome=o, model='KRAS-mut vs WTall')
T = pd.DataFrame(rows)
# ---- H5 within EGFR-mut: DUSP6 vs pERK and ETV
R = []
E = A[A.group == 'EGFR']
def pcor(d, x, y, cov=('purity',)):
    d = d[[x, y, *cov]].dropna()
    if len(d) < 10: return np.nan, np.nan, len(d)
    rx = sm.OLS(stats.rankdata(d[x]), sm.add_constant(d[list(cov)])).fit().resid
    ry = sm.OLS(stats.rankdata(d[y]), sm.add_constant(d[list(cov)])).fit().resid
    r, p = stats.pearsonr(rx, ry); return r, p, len(d)
for sub, D in [('EGFR-mut', E), ('RAS-path', A[A.group == 'RAS']), ('WTall', A[A.group == 'WT']), ('all LUAD', A)]:
    for x in ['rna_DUSP6', 'prot_DUSP6', 'DUSP6_protBroad', 'rna_DUSP4', 'prot_DUSP4', 'rna_SPRY4', 'rna_SPRED1', 'stage_RTK']:
        for y in ['stage_ERK', 'pERK_broad', 'stage_ERKSUB', 'rna_ETV4', 'rna_ETV5', 'ETVc']:
            if x not in A or A[x].notna().sum() < 10: continue
            r, p, n = pcor(D, x, y); R.append(dict(subset=sub, x=x, y=y, rho_partial_purity=r, p=p, n=n))
    for y in ['rna_ETV4', 'rna_ETV5', 'ETVc', 'stage_ERKSUB']:
        for x in ['stage_ERK', 'pERK_broad']:
            r, p, n = pcor(D, x, y); R.append(dict(subset=sub, x=x, y=y, rho_partial_purity=r, p=p, n=n))
Rt = pd.DataFrame(R)
T.to_csv(OUT + 'tables/T2_robustness_panel.tsv', sep='\t', index=False); Rt.to_csv(OUT + 'tables/T3_within_group_correlations.tsv', sep='\t', index=False)
# ---- ETV/pERK coupling slope by group (does EGFR-mut have more ETV per unit pERK?)
S = []
for x in ['stage_ERK', 'pERK_broad']:
    d = A[A.group.isin(['EGFR', 'RAS', 'WT'])][['ETVc', x, 'group', 'purity', 'prolif']].dropna()
    X = pd.DataFrame({'pERK': z(d[x]), 'EGFR': (d.group == 'EGFR') * 1., 'RAS': (d.group == 'RAS') * 1., 'purity': z(d.purity), 'prolif': z(d.prolif)})
    X['EGFRxpERK'] = X.EGFR * X.pERK; X['RASxpERK'] = X.RAS * X.pERK
    r = sm.OLS(z(d.ETVc).values, sm.add_constant(X)).fit(cov_type='HC3')
    for k in ['pERK', 'EGFR', 'RAS', 'EGFRxpERK', 'RASxpERK']:
        S.append(dict(x=x, term=k, beta=r.params[k], lo=r.conf_int().loc[k, 0], hi=r.conf_int().loc[k, 1], p=r.pvalues[k], n=len(d)))
S = pd.DataFrame(S); S.to_csv(OUT + 'tables/T4_ETV_on_pERK_interaction.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
c = ['test', 'outcome', 'model', 'contrast', 'beta', 'lo', 'hi', 'p', 'n_g', 'n_ref']
print(T[(T.test == 'panel') & (T.model == 'base (vs WTall)')][c].round(3).to_string())
print(T[(T.test == 'panel') & (T.model == 'EGFR vs RAS-path')][['outcome', 'beta', 'lo', 'hi', 'p']].round(3).to_string())
print(T[T.test != 'panel'][c].round(3).to_string())
print(Rt[Rt.subset == 'EGFR-mut'].round(3).to_string())
print(S.round(3).to_string())
