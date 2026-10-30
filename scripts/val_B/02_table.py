import json, pandas as pd, numpy as np, sys
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
rows = []
for p in ['CPTAC-3', 'CPTAC-2']:
    for c in json.load(open(f'/workspace/cache/val_B/dl/gdc_{p}_cases.json')):
        for s in c.get('samples', []):
            if s.get('tissue_type') != 'Tumor': continue
            rows.append(dict(case_id=c['submitter_id'], proj=p, excision_min=s.get('time_between_excision_and_freezing'), clamping_min=s.get('time_between_clamping_and_freezing')))
d = pd.DataFrame(rows).astype({'excision_min': float, 'clamping_min': float})
g = d.groupby('case_id').agg(excision_min=('excision_min', 'median'), clamping_min=('clamping_min', 'median'), proj=('proj', 'first'))
g.to_csv('/workspace/results/val_B/gdc_ischaemia_all.tsv', sep='\t')
s = cl.samples(); s = s[(s.sample_type == 'Tumor') & s.has_phospho]
cases = s.drop_duplicates('case_id')[['case_id', 'cohort']]
m = cases.merge(g, on='case_id', how='left')
cov = m.groupby('cohort').agg(n_tumour_cases=('case_id', 'size'), in_gdc=('proj', lambda x: x.notna().sum()), n_excision=('excision_min', lambda x: x.notna().sum()),
                              n_clamp=('clamping_min', lambda x: x.notna().sum()), med_excision=('excision_min', 'median'), q1=('excision_min', lambda x: x.quantile(.25)), q3=('excision_min', lambda x: x.quantile(.75)),
                              min_ex=('excision_min', 'min'), max_ex=('excision_min', 'max'), n_unique=('excision_min', 'nunique'))
cov.to_csv('/workspace/results/val_B/T1_coverage.tsv', sep='\t'); print(cov.round(1).to_string())
