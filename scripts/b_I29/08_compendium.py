"""I29 H8: compendium - ERK activation-loop sites under EGFR-TKI in EGFR-mutant lines (PC9 decryptM dose series; PMID 25404012 lines)."""
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, '/workspace/scripts/infra_compendium'); import load
OUT = '/workspace/results/b_I29/'
C = load.conditions()
m = C.cell_line_norm.fillna('').str.upper()
sel = C[m.isin(['PC9', 'NCIH1975', 'NCIH3255', 'HCC827', 'HCC4006', 'A549', 'NCIH358', 'NCIH1299'])]
print(sel[sel.study_pmid.isin(['25404012'])][['conditionID', 'cell_line_norm', 'perturbation_class', 'agents', 'stimulus', 'time_min', 'dose_nM', 'contrast_reversed']].to_string())
SITES = {'MAPK3_T202': 'P27361_T202', 'MAPK3_Y204': 'P27361_Y204', 'MAPK1_T185': 'P28482_T185', 'MAPK1_Y187': 'P28482_Y187', 'EGFR_Y1197': 'P00533_Y1197', 'EGFR_Y1092': 'P00533_Y1092', 'EGFR_Y1172': 'P00533_Y1172',
         'RSK1_S380': 'Q15418_S380', 'RSK3_S386': 'P51812_S386', 'ETV? none': 'none'}
X = load.sites_matrix(conditions=list(sel.conditionID), sites=[v for v in SITES.values() if v != 'none'])
psp = pd.read_csv('/workspace/data/Kinase_Substrate_Dataset', sep='\t', skiprows=3, encoding='latin-1')
psp = psp[(psp.KIN_ORGANISM == 'human') & (psp.SUB_ORGANISM == 'human') & psp.KIN_ACC_ID.isin(['P28482', 'P27361'])]
es = sorted(set(psp.SUB_ACC_ID + '_' + psp.SUB_MOD_RSD) - {v for v in SITES.values()})
Y = load.sites_matrix(conditions=list(sel.conditionID), sites=es)
R = X.T.rename(columns={v: k for k, v in SITES.items()})
R['ERKSUB_mean'] = Y.mean(0); R['ERKSUB_n'] = Y.notna().sum(0)
R = sel.set_index('conditionID')[['study_pmid', 'cell_line_norm', 'agents', 'stimulus', 'time_min', 'dose_nM', 'contrast_reversed']].join(R)
R.to_csv(OUT + 'tables/T10_compendium_conditions.tsv', sep='\t')
rows = []
for (cl_, ag), g in R[R.agents.notna()].groupby(['cell_line_norm', 'agents']):
    g = g.sort_values('dose_nM')
    for s in ['MAPK3_T202', 'MAPK3_Y204', 'MAPK1_T185', 'MAPK1_Y187', 'EGFR_Y1197', 'EGFR_Y1092', 'RSK1_S380', 'ERKSUB_mean']:
        d = g[['dose_nM', s]].dropna()
        if len(d) >= 4:
            rho, p = stats.spearmanr(d.dose_nM, d[s])
            rows.append(dict(cell_line=cl_, agents=ag, site=s, n=len(d), rho_dose=rho, p=p, log2fc_top2dose=d[s].iloc[-2:].mean(), log2fc_low2dose=d[s].iloc[:2].mean()))
T = pd.DataFrame(rows); T.to_csv(OUT + 'tables/T11_compendium_tki_dose.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(T.round(2).to_string())
print(R[R.study_pmid == '25404012'].round(2).to_string())
