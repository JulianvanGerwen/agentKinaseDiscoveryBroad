"""I29: discovery LUAD tumours re-measured inside the Confirmatory plexes (CDAP): does the CDAP measurement reproduce the discovery EGFR effect in the SAME tumours?"""
import re, json, numpy as np, pandas as pd
exec(open('/workspace/scripts/b_I29/07_replication_cptac_conf.py').read().split("disc = set(")[0])
disc = pd.read_parquet('/workspace/results/b_I27/tables/alleles.parquet').query('cohort=="LUAD"').case_id
meta = pd.DataFrame({k: dict(case=v[0], tissue=v[1]) for k, v in al.items()}).T
keep = [i for i, c in enumerate(PM.columns) if c in meta.index and meta.loc[c, 'tissue'] == 'Tumor' and meta.loc[c, 'case'] in set(disc)]
T = PM.iloc[:, keep].T; T.index = [meta.loc[c, 'case'] for c in T.index]; T = T.groupby(level=0).mean()
print('discovery tumours in confirmatory plexes:', len(T))
import sys; sys.path.insert(0, '/workspace/scripts/b_I29'); from lib import load
A = load('LUAD'); A = A.join(T[[c for c in T.columns if c.startswith('MAPK')]].add_prefix('conf_'))
for c in ['conf_MAPK3:y204', 'conf_MAPK1:y187', 'conf_MAPK3:t202y204']:
    d = A[[c, 'MAPK3_pY_raw', 'stage_ERK', 'group']].dropna(subset=[c])
    print(c, len(d), d.group.value_counts().to_dict(), 'r vs discovery Broad MAPK3_pY_raw', round(d[[c, 'MAPK3_pY_raw']].corr().iloc[0, 1], 2), 'r vs stage_ERK', round(d[[c, 'stage_ERK']].corr().iloc[0, 1], 2))
    print('  group medians conf:', d.groupby('group')[c].median().round(2).to_dict(), ' discovery Broad pY_raw:', d.groupby('group').MAPK3_pY_raw.median().round(2).to_dict())
