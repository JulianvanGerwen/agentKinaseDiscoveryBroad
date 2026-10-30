"""val_E step 1: per-tumour cold-ischaemia score (audit_I10 definition) for all CPTAC cohorts; site-ID proxy; circularity check."""
import importlib.util, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
sp = importlib.util.spec_from_file_location('c10', '/workspace/scripts/audit_I10/common.py'); c10 = importlib.util.module_from_spec(sp); sys.path.insert(0, '/workspace/scripts/audit_I10'); sp.loader.exec_module(c10)
OUT = '/workspace/results/val_E/'
ISCH = c10.ischaemia_set(200)
PATH = {'P00533','P04626','P08581','P29353','Q13480','Q06124','P27361','P28482','P21860'}
pathsites = [s for s in ISCH.index if s.split('_')[0] in PATH]
print('ISCH sites on pathway proteins (EGFR/ERBB2/MET/SHC1/GAB1/PTPN11/ERK1/2):', pathsites)
mapk_like = [s for s in ISCH.index if s.split('_')[0] in {'P27361','P28482','P15056','Q02750','P36507','P49137','Q15418','P51812','P04792'}]
print('MAPK-axis sites in ISCH:', mapk_like)
rows = []
for variant in ['input_tc_padj', 'input_tc']:
    for coh in c10.COH:
        Z = c10.standardise(c10.cohort(coh, variant=variant))
        s = c10.score(Z, ISCH)
        s2 = c10.score(Z, ISCH[~ISCH.index.str.split('_').str[0].isin(PATH)])
        d = pd.DataFrame({'isch': s, 'isch_noPath': s2}); d['cohort'] = coh; d['variant'] = variant
        d['n_isch_sites'] = Z.reindex(ISCH.index.intersection(Z.index)).notna().sum(0)
        rows.append(d)
A = pd.concat(rows); A['case_id'] = [i[:-2] for i in A.index]; A['site_prefix'] = A.case_id.str.extract(r'^(\d+)')[0]
A.to_parquet('/workspace/cache/val_E/isch_scores.parquet'); A.to_csv(OUT + 'ischaemia_scores.tsv', sep='\t')
print(A.groupby(['variant','cohort']).agg(n=('isch','count'), sd=('isch','std'), r_noPath=('isch', lambda x: np.corrcoef(x.dropna(), A.loc[x.dropna().index.to_series().index, 'isch_noPath'].reindex(x.dropna().index).fillna(0))[0,1] if False else np.nan)).drop(columns='r_noPath').round(3))
print(A.groupby(['variant','cohort']).apply(lambda g: g[['isch','isch_noPath']].corr().iloc[0,1]).round(3).unstack(0))
print(A[A.variant=='input_tc_padj'].groupby('cohort').site_prefix.agg(lambda x: x.nunique()))
