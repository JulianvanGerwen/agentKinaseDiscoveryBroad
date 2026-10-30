"""Fetch GDC sample-level time_between_excision_and_freezing / clamping for CPTAC-2 and CPTAC-3 (all cases)."""
import json, requests, sys
DL = '/workspace/cache/val_B/dl/'
for proj in ['CPTAC-3', 'CPTAC-2']:
    out = []; frm = 0
    while True:
        r = requests.post('https://api.gdc.cancer.gov/cases', timeout=120, json=dict(
            filters={"op": "in", "content": {"field": "project.project_id", "value": [proj]}},
            fields='submitter_id,project.project_id,primary_site,disease_type,samples.submitter_id,samples.sample_type,samples.tissue_type,samples.time_between_excision_and_freezing,samples.time_between_clamping_and_freezing,samples.preservation_method,samples.method_of_sample_procurement,samples.tumor_descriptor',
            size=200, **{'from': frm}))
        r.raise_for_status(); d = r.json()['data']; out += d['hits']; frm += 200
        if frm >= d['pagination']['total']: break
    json.dump(out, open(DL + f'gdc_{proj}_cases.json', 'w')); print(proj, len(out))
