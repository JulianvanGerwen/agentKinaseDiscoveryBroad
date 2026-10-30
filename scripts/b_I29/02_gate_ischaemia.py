"""I29 Gate G: cold-ischaemia / collection-site confound for EGFR-mut pERK-low / ETV-high dissociation (LUAD)."""
import numpy as np, pandas as pd
from scipy import stats
from lib import *
A = load('LUAD')
pd.set_option('display.width', 250)
# ---- G1: ischaemia proxies by group
rows = []
for v in ['ISCH', 'ISCH_noERK', 'isch_valE', 'isch_noPath_valE', 'excision_min', 'clamping_min', 'log_excision', 'rna_FOS', 'rna_EGR1', 'rna_JUN']:
    r = ols(A, v, ['purity', 'prolif']); rows += row(r, test='G1 group difference', outcome=v, model='base')
    kw = stats.kruskal(*[A.loc[A.group == g, v].dropna() for g in ['EGFR', 'RAS', 'WT']])
    rows.append(dict(test='G1 Kruskal', outcome=v, model='KW 3 groups', p=kw.pvalue,
                     beta=np.nan, contrast='; '.join(f"{g} median {A.loc[A.group == g, v].median():.3f}" for g in ['EGFR', 'RAS', 'WT'])))
# country x group
ct = pd.crosstab(A.country, A.group); print(ct)
rows.append(dict(test='G1 country x group', outcome='country', model='chi2', p=stats.chi2_contingency(ct)[1], contrast=ct.to_dict().__str__()))
# ---- validity of ISCH as ischaemia proxy: correlation with real cold ischaemia time
for v in ['ISCH', 'ISCH_noERK', 'isch_valE', 'stage_ERK', 'pERK_broad', 'stage_ERKSUB', 'rna_FOS', 'rna_EGR1', 'ETVc']:
    d = A[[v, 'excision_min']].dropna(); rho, p = stats.spearmanr(d[v], d.excision_min)
    rows.append(dict(test='G0 Spearman vs excision->freezing time', outcome=v, model='all LUAD', beta=rho, p=p, n_g=len(d)))
    d = A[[v, 'clamping_min']].dropna(); rho, p = stats.spearmanr(d[v], d.clamping_min)
    rows.append(dict(test='G0 Spearman vs clamping->freezing time', outcome=v, model='all LUAD', beta=rho, p=p, n_g=len(d)))
for v in ['stage_ERK', 'pERK_broad', 'stage_ERKSUB', 'ETVc']:
    d = A[[v, 'ISCH_noERK']].dropna(); rho, p = stats.spearmanr(d[v], d.ISCH_noERK)
    rows.append(dict(test='G0 Spearman vs ISCH_noERK', outcome=v, model='all LUAD', beta=rho, p=p, n_g=len(d)))
# ---- G2: contrasts under adjustment
B = ['purity', 'prolif']
MODELS = {'base': (B, ()), '+ISCH_noERK': (B + ['ISCH_noERK'], ()), '+ISCH(all 200)': (B + ['ISCH'], ()), '+isch_noPath(valE)': (B + ['isch_noPath_valE'], ()),
          '+log excision time': (B + ['log_excision_f'], ()), '+country': (B, ('country',)),
          'FULL: +ISCH_noERK+country+log excision': (B + ['ISCH_noERK', 'log_excision_f'], ('country',)),
          'FULL + TMT plex FE': (B + ['ISCH_noERK', 'log_excision_f'], ('country', 'EXPERIMENT'))}
OUTS = ['stage_ERK', 'pERK_broad', 'MAPK1_pTpY', 'MAPK3_pTpY', 'stage_ERKSUB', 'stage_RSK', 'MAPKout_fp', 'stage_RTK', 'ETVc', 'rna_ETV4', 'rna_ETV5', 'rna_SPRED1', 'rna_DUSP6', 'rna_DUSP4']
for o in OUTS:
    for m, (cv, fe) in MODELS.items():
        nperm = 2000 if (o in ('stage_ERK', 'ETVc') and m in ('base', 'FULL: +ISCH_noERK+country+log excision')) else 0
        rows += row(ols(A, o, cv, fe=fe, nperm=nperm), test='G2 adjusted contrast', outcome=o, model=m)
    for nm, msk in [('within Vietnam', A.country == 'Vietnam'), ('within non-Vietnam', A.country != 'Vietnam'), ('within Asian (VN+CN)', A.country.isin(['Vietnam', 'China'])), ('within USA', A.country == 'Usa'),
                    ('excision time known, <= median', A.excision_min <= A.excision_min.median()), ('excision time > median', A.excision_min > A.excision_min.median())]:
        rows += row(ols(A, o, B + ['ISCH_noERK'], mask=msk), test='G2 stratum', outcome=o, model=nm)
T = pd.DataFrame(rows); T.to_csv(OUT + 'tables/T1_gate_ischaemia.tsv', sep='\t', index=False)
show = T[T.test.str.startswith(('G0', 'G1'))]
print(show[['test', 'outcome', 'model', 'contrast', 'beta', 'p', 'n_g']].round(3).to_string())
g = T[(T.test.str.startswith('G2')) & (T.contrast == 'EGFR_vs_ref')]
print(g.pivot_table(index='outcome', columns='model', values='beta', sort=False).round(2).to_string())
print(g.pivot_table(index='outcome', columns='model', values='p', sort=False).round(3).to_string())
print(T[T.p_perm.notna()][['outcome', 'model', 'contrast', 'beta', 'lo', 'hi', 'p', 'p_perm']].round(3).to_string()) if 'p_perm' in T else None
