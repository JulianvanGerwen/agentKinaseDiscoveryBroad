"""I27 step 4: follow-ups. (a) collapse degenerate MAP2K family scores; (b) smoking confound in LUAD;
(c) direct MAPK/JNK/p38 activation-loop sites; (d) P1 RAF-footprint follow-up in PDAC."""
import numpy as np, pandas as pd, common as C

OUT = '/workspace/results/b_I27/tables/'
A = pd.read_parquet(OUT + 'alleles.parquet')
X = pd.read_parquet(OUT + 'outcomes_tc_padj.parquet')
D = X.merge(A, on=['case_id', 'cohort'], how='left')
D['kin_MAP2Kfam'] = D[['kin_MAP2K3', 'kin_MAP2K4', 'kin_MAP2K5', 'kin_MAP2K6', 'kin_MAP2K7']].mean(1)
D['kin_MEK12'] = D[['kin_MAP2K1', 'kin_MAP2K2']].mean(1)

# direct activation-loop sites (MAPK1/3 TEY, MAPK8/9 TPY, MAPK14 TGY, MAP2K1/2 S218/222/S222/226)
LOOP = {'ERK_TEY': ['P28482_Y187', 'P28482_T185', 'P27361_T202', 'P27361_Y204'],
        'JNK_TPY': ['P45983_T183', 'P45983_Y185', 'P45984_T183', 'P45984_Y185'],
        'p38_TGY': ['Q16539_T180', 'Q16539_Y182'],
        'MEK_loop': ['Q02750_S218', 'Q02750_S222', 'P36507_S222', 'P36507_S226']}
sit = []
for coh in C.COHORTS:
    ph = C.cl.matrix(coh, 'input_tc_padj')
    tt = [s for s in C.cl.tumour_ids(coh) if s in ph.columns]
    df = pd.DataFrame(index=tt)
    for k, ss in LOOP.items():
        ss2 = [s for s in ss if s in ph.index]
        df['loop_' + k] = ph.loc[ss2, tt].mean(0) if ss2 else np.nan
    df['cohort'] = coh
    sit.append(df)
S = pd.concat(sit)
lc = [c for c in S.columns if c.startswith('loop_')]
S[lc] = S.groupby('cohort')[lc].transform(lambda s: (s - s.mean()) / s.std(ddof=1))
S['case_id'] = [s[:-2] for s in S.index]
D = D.merge(S[['case_id'] + lc], on='case_id', how='left')
print('loop site coverage:', D[lc].notna().mean().round(2).to_dict())

# smoking covariate for LUAD
sm = C.cl.table('LUAD', 'subtypes_cbioportal')
col = 'SMOKING_SCORE_WGS' if 'SMOKING_SCORE_WGS' in sm.columns else 'PACK_YEARS_SMOKED'
smk = pd.DataFrame({'case_id': sm.case_id, 'smoke': pd.to_numeric(sm[col], errors='coerce'),
                    'smoke_status': sm['SMOKING_STATUS'].astype(str) if 'SMOKING_STATUS' in sm.columns else 'NA'})
D = D.merge(smk.drop_duplicates('case_id'), on='case_id', how='left')
D.loc[D.cohort == 'LUAD', 'smoke'] = D.loc[D.cohort == 'LUAD', 'smoke'].fillna(D.loc[D.cohort == 'LUAD', 'smoke'].median())

rows = []
L = D[D.cohort == 'LUAD'].copy()
L['EGFRmut'] = np.where(L.EGFR_allele.isin(['ex19del', 'L858R']), 'EGFRmut', np.where(L.EGFR_allele == 'WT', 'WT', 'x'))
Y = ['kin_MAP2Kfam', 'kin_MEK12', 'kin_MAPK1', 'loop_ERK_TEY', 'loop_JNK_TPY', 'loop_p38_TGY', 'loop_MEK_loop',
     'stage_ERK', 'stage_ERKSUB', 'stage_AKTSUB', 'stage_RTK', 'stage_ADAPT', 'fp_KI_MEK', 'fp_MAPKaxis']
for a, b, gcol in [('G12C', 'WT', 'KRAS_allele'), ('L858R', 'WT', 'EGFR_allele'),
                   ('ex19del', 'WT', 'EGFR_allele'), ('ex19del', 'L858R', 'EGFR_allele'),
                   ('EGFRmut', 'WT', 'EGFRmut')]:
    for y in Y:
        for lab, extra in [('base', ()), ('smoke_adj', ('smoke',))]:
            r = C.contrast(L, y, gcol, a, b, extra=extra, nperm=2000, seed=3)
            r.update(model=lab, cohort='LUAD'); rows.append(r)
# PDAC P1 follow-up: RAF / MEK / PI3K footprints and MAPK loop
P = D[D.cohort == 'PDAC'].copy()
P['K_P1'] = np.where(P.KRAS_allele == 'G12R', 'G12R', np.where(P.KRAS_allele.isin(['G12D', 'G12V']), 'G12DV', 'x'))
for y in ['fp_KI_RAF', 'fp_KI_MEK', 'fp_KI_PI3K', 'fp_PI3Kaxis', 'fp_MAPKaxis', 'loop_ERK_TEY', 'loop_MEK_loop',
          'stage_ERK', 'stage_ERKSUB', 'stage_AKTSUB', 'kin_MAP2Kfam']:
    for a, b in [('G12R', 'G12DV'), ('G12R', 'G12D'), ('G12V', 'G12D')]:
        g = 'K_P1' if b == 'G12DV' else 'KRAS_allele'
        r = C.contrast(P, y, g, a, b, nperm=2000, seed=4)
        r.update(model='PDAC_followup', cohort='PDAC'); rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(OUT + 'T6_followup.tsv', sep='\t', index=False)
print(R[R.p < 0.05][['cohort', 'a', 'b', 'y', 'model', 'n_a', 'n_b', 'beta', 'lo', 'hi', 'p']].to_string())
