"""Per time series: site-level peak, direction and onset.
Series kept: perturbations with a defined PRIMARY direction (inhibitors: down; activating stimuli: up).
Time points <=0.5 min dropped, replicate time points averaged. Sites must be measured at every time point of the series.
Onset (per threshold thr): first time index with dir*x > thr, where dir = sign of the peak (max |x|).
Also half-max onset: first index with dir*x >= 0.5*|peak| (magnitude-free).
Outputs: cache/series_meta.tsv, cache/onset_long.parquet"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_compendium'); import load
OUT = '/workspace/results/b_I05/cache/'
C = load.conditions()
T = C[C.time_series_id.notna() & ~C.qc_flag.isin(['scale_suspect', 'duplicate_vector'])].copy()

# ---- series annotation: primary direction, replication class (audit_I19 classes), target / expected kinases
UPSUB = {'RTK_ligand', 'serum', 'TGFb', 'TNF', 'interferon', 'interleukin', 'IFNG_LPS', 'GPCR_ligand', 'immune_receptor',
         'PKC_activator', 'phosphatase_inhibitor', 'DNA_damage', 'UV', 'antibody_CD20', 'receptor_agonist'}
INH_CLASS = {'vemurafenib': 'RAF', 'gefitinib': 'EGFR', 'PF04217903': 'MET_ALK', 'MK2206': 'AKT', 'rapamycin': 'rapalog', 'ATRi': 'DDR',
             'ATMi': 'DDR', 'DNAPKi': 'DDR', 'THZ1': 'CDK', 'CDK12IN3': 'CDK', 'dinaciclib': 'CDK', 'tofacitinib': 'JAK', 'gefitinib;vemurafenib': 'EGFR'}
# pathway effectors added to the inhibitor's primary targets (targets that act mainly via one downstream kinase)
EFFECT = {'BRAF': ['MAPK1', 'MAPK3'], 'EGFR': ['MAPK1', 'MAPK3'], 'MET': ['MAPK1', 'MAPK3'], 'MTOR': ['RPS6KB1'], 'ATR': ['ATM', 'PRKDC'],
          'ATM': ['ATR', 'PRKDC'], 'PRKDC': ['ATM', 'ATR']}
STIM_K = {'LPA': 'MAPK1;MAPK3;PRKCA', 'thrombin': 'MAPK1;MAPK3;PRKCA', 'endothelin': 'MAPK1;MAPK3;PRKCA', 'angiotensin II': 'MAPK1;MAPK3;PRKCA',
          'serum': 'MAPK1;MAPK3;RPS6KA1;AKT1', 'NCS': 'ATM;ATR;PRKDC', 'GLP1': 'PRKACA', 'TGF-beta': 'TGFBR1;MAPK1;MAPK3'}
GPCR = {'LPA', 'thrombin', 'endothelin', 'angiotensin II', 'PGE2', 'ADP'}
rows = []
for sid, d in T.groupby('time_series_id'):
    r = d.iloc[0]; sub = r.perturbation_subclass; prim = None; rcls = None; kin = ''
    if r.perturbation_class == 'kinase_inhibitor' and isinstance(r.kinase_inhibitors, str):
        if sid == 'TS0042': continue  # tofacitinib cocult vs monocult: confounded with co-culture
        if r.cell_line_norm == 'HCT116' and r.kinase_inhibitors == 'vemurafenib': prim = 'para'  # BRAF-WT/KRAS: paradoxical activation
        else: prim = 'down'
        rcls = INH_CLASS.get(r.kinase_inhibitors)
        tg = str(r.inhibitor_primary_targets).split(';')
        kin = ';'.join(dict.fromkeys(tg + sum([EFFECT.get(t, []) for t in tg], [])))
    elif r.perturbation_class == 'genetic' or r.perturbation_class == 'time_course':
        continue
    elif sub in UPSUB or r.perturbation_class in ('growth_factor', 'cytokine'):
        if sid in ('TS0014', 'TS0015', 'TS0032', 'TS0033', 'TS0072', 'TS0074'): continue  # platelet <=1 min / unclassified FGFR ligands
        prim = 'up'
        st = r.stimulus if isinstance(r.stimulus, str) else ''
        if sub == 'DNA_damage': st = 'NCS'
        if sub == 'receptor_agonist' and 'glutide' in str(r.condition_raw) or 'Danuglipron' in str(r.condition_raw): st = 'GLP1'
        kin = r.stimulus_expected_up_kinases if isinstance(r.stimulus_expected_up_kinases, str) else STIM_K.get(st, '')
        rcls = ('GF_RTK' if r.perturbation_class == 'growth_factor' else 'TNF' if st == 'TNF' else 'GPCR' if st in GPCR else
                'DNA_damage' if sub in ('DNA_damage', 'UV') else 'PMA' if st == 'PMA' else 'PPase_inh' if sub == 'phosphatase_inhibitor'
                else 'cytokine_other' if r.perturbation_class == 'cytokine' else None)
    elif sub == 'starvation' and sid == 'TS0013':
        prim = 'down'; kin = 'MTOR;RPS6KB1'
    else:
        continue
    rows.append(dict(series=sid, study=r.study_pmid, cell=r.cell_line_norm, pclass=r.perturbation_class, sub=sub, label=r.condition_raw,
                     agent=r.kinase_inhibitors if isinstance(r.kinase_inhibitors, str) else (r.stimulus if isinstance(r.stimulus, str) else r.agents),
                     prim=prim, rcls=rcls, kinases=kin))
SM = pd.DataFrame(rows)
T = T[T.time_series_id.isin(SM.series) & (T.time_min > 0.5)]
X = load.sites_matrix(T.conditionID.tolist())
X = X * np.where(T.set_index('conditionID').loc[X.columns, 'contrast_reversed'].astype(bool), -1.0, 1.0)
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'residue', 'psp_kinases_human'], low_memory=False).set_index('site')
MP = pd.read_parquet(OUT + 'motif_pct.parquet')
THR = [0.5, 1.0]
long = []; meta = []; VV = {}
for s in SM.itertuples():
    d = T[T.time_series_id == s.series]
    times = np.sort(d.time_min.unique())
    if len(times) < 2: continue
    V = np.column_stack([X[d.conditionID[d.time_min == t]].mean(axis=1).values for t in times])
    ok = ~np.isnan(V).any(1)
    V = V[ok]; sites = X.index[ok]
    if len(sites) < 300: continue
    pk_i = np.abs(V).argmax(1); peak = V[np.arange(len(V)), pk_i]; dirn = np.sign(peak)
    dv = V * dirn[:, None]
    rec = dict(series=s.series, site=sites, peak=peak.astype(np.float32), peak_idx=pk_i.astype(np.int8), last=V[:, -1].astype(np.float32), first=V[:, 0].astype(np.float32))
    for thr in THR:
        hit = dv > thr
        rec[f'on{thr}'] = np.where(hit.any(1), hit.argmax(1), -1).astype(np.int8)
    hm = dv >= 0.5 * np.abs(peak)[:, None]
    rec['on_half'] = hm.argmax(1).astype(np.int8)
    L = pd.DataFrame(rec)
    L['dir'] = np.where(L.peak > 0, 'up', 'down')
    L['primary'] = (L.dir == s.prim) if s.prim in ('up', 'down') else np.nan
    ks = [k for k in s.kinases.split(';') if k in MP.columns]
    L['res'] = S.residue.reindex(L.site).values
    L['tgt_pct'] = MP.reindex(L.site)[ks].max(axis=1).values if ks else np.nan
    psp = S.psp_kinases_human.reindex(L.site).fillna('').values
    L['tgt_psp'] = [bool(set(p.split(';')) & set(ks)) for p in psp] if ks else False
    long.append(L); VV[s.series] = (np.array(sites), V.astype(np.float32), times)
    meta.append(dict(s._asdict(), n_tp=len(times), times=','.join(f'{t:g}' for t in times), n_sites=len(L), kin_used=';'.join(ks),
                     n_resp1=int((np.abs(peak) > 1).sum())))
SM2 = pd.DataFrame(meta).drop(columns='Index')
SM2.to_csv(OUT + 'series_meta.tsv', sep='\t', index=False)
LL = pd.concat(long, ignore_index=True)
LL.to_parquet(OUT + 'onset_long.parquet')
import pickle; pickle.dump(VV, open(OUT + 'series_V.pkl', 'wb'))
print(SM2.to_string()); print(LL.shape)
