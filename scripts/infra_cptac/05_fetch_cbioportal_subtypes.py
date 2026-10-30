#!/usr/bin/env python3
"""Fetch per-sample/patient clinical attributes (incl. molecular subtypes) of the CPTAC cohort-paper studies on
cBioPortal (public API) and write a wide table per cohort: /workspace/results/infra_cptac/cohorts/<C>/subtypes_cbioportal.parquet
(case_id + all attributes as strings). Raw JSON cached in /workspace/cache/cptac/cbioportal/."""
import json, os, urllib.request
import pandas as pd

STUDIES = {'BRCA': 'brca_cptac_2020', 'COAD': 'coad_cptac_2019', 'UCEC': 'ucec_cptac_2020', 'GBM': 'gbm_cptac_2021',
           'LUAD': 'luad_cptac_2020', 'LSCC': 'lusc_cptac_2021', 'PDAC': 'paad_cptac_2021'}
CACHE = '/workspace/cache/cptac/cbioportal'; os.makedirs(CACHE, exist_ok=True)
OUT = '/workspace/results/infra_cptac/cohorts'
API = 'https://www.cbioportal.org/api'


def get(url, fn):
    p = f'{CACHE}/{fn}'
    if not os.path.exists(p):
        with urllib.request.urlopen(url, timeout=120) as r:
            open(p, 'wb').write(r.read())
    return json.load(open(p))


summary = {}
for c, sid in STUDIES.items():
    frames = []
    for typ in ('PATIENT', 'SAMPLE'):
        d = get(f'{API}/studies/{sid}/clinical-data?clinicalDataType={typ}&projection=SUMMARY&pageSize=100000', f'{sid}_{typ}.json')
        df = pd.DataFrame(d)
        if df.empty: continue
        key = 'patientId'
        w = df.pivot_table(index=key, columns='clinicalAttributeId', values='value', aggfunc='first')
        frames.append(w)
    w = pd.concat(frames, axis=1)
    w = w.loc[:, ~w.columns.duplicated()]
    w.index = w.index.astype(str).str.replace(r'^X(?=\d)', '', regex=True); w.index.name = 'case_id'   # R-mangled BRCA IDs (X01BR001)
    w.insert(0, 'cbioportal_study', sid)
    w.reset_index().astype(str).to_parquet(f'{OUT}/{c}/subtypes_cbioportal.parquet', index=False)
    ours = pd.read_parquet(f'{OUT}/{c}/samples.parquet').case_id.unique()
    summary[c] = {'study': sid, 'n_cases': len(w), 'n_overlap_with_phospho_cases': int(pd.Index(ours).isin(w.index).sum()),
                  'columns': list(w.columns)}
    print(c, sid, len(w), summary[c]['n_overlap_with_phospho_cases'], [x for x in w.columns if any(k in x.upper() for k in ('SUBTYPE', 'PAM50', 'CLUSTER', 'NMF', 'MSI', 'IMMUNE', 'CLASS'))])
json.dump(summary, open('/workspace/results/infra_cptac/subtypes_cbioportal_summary.json', 'w'), indent=1)
