"""I27 step 2: pre-registered tests P1-P5."""
import numpy as np, pandas as pd, common as C

OUT = '/workspace/results/b_I27/tables/'
A = pd.read_parquet(OUT + 'alleles.parquet')
res = []


def prep(variant='tc_padj'):
    X = pd.read_parquet(OUT + f'outcomes_{variant}.parquet')
    D = X.merge(A, on=['case_id', 'cohort'], how='left')
    # subtype covariates
    sub = []
    for coh in C.COHORTS:
        try:
            s = C.cl.table(coh, 'subtypes_cbioportal')
        except Exception:
            continue
        col = None
        for cand in ['PAM50', 'NMF_CLUSTER', 'MSI_STATUS', 'HISTOLOGIC_TYPE', 'CANCER_TYPE_DETAILED']:
            if cand in s.columns:
                col = cand; break
        if col is None: continue
        sub.append(pd.DataFrame({'case_id': s.case_id, 'subtype': s[col].astype(str)}))
    if sub:
        S = pd.concat(sub).drop_duplicates('case_id')
        D = D.merge(S, on='case_id', how='left')
        D['subtype'] = D['subtype'].fillna('NA')
    return D


D = prep('tc_padj')
Dtc = prep('tc')
STAGES = ['stage_RTK', 'stage_EGFRpY', 'stage_ADAPT', 'stage_ERK', 'stage_RSK', 'stage_ERKSUB',
          'stage_AKTSUB', 'stage_MTORC1', 'stage_S6']
AKTK = ['kin_AKT1', 'kin_AKT2', 'kin_MTOR', 'kin_RPS6KB1', 'kin_SGK1', 'kin_SGK3', 'kin_PDPK1',
        'kin_MAPK1', 'kin_MAPK3', 'kin_RPS6KA1', 'kin_RPS6KA3']
SITES = ['site_PRAS40_T246', 'site_NDRG1_SGK', 'site_GSK3B_S9', 'site_S6_S235']
FPS = ['fp_KI_PI3K', 'fp_KI_AKT', 'fp_KI_MEK', 'fp_KI_mTORcat', 'fp_KI_EGFR', 'fp_EGF', 'fp_GF_all']
ALLY = STAGES + AKTK + SITES + FPS


def run(tag, df, cohorts, gcol, a, b, ycols, extra=(), note=''):
    rows = []
    for coh in cohorts:
        d = df[df.cohort == coh]
        for y in ycols:
            r = C.contrast(d, y, gcol, a, b, extra=extra)
            r.update(pred=tag, cohort=coh, note=note)
            rows.append(r)
    res.extend(rows)
    return pd.DataFrame(rows)


# ---------- P1: KRAS G12R vs G12D+G12V in PDAC ----------
d = D.copy()
d['K_P1'] = np.where(d.KRAS_allele == 'G12R', 'G12R',
                     np.where(d.KRAS_allele.isin(['G12D', 'G12V']), 'G12DV', 'other'))
# primary: AKTSUB adjusted for ERK ("matched MAPK output")
r = C.contrast(d[d.cohort == 'PDAC'], 'stage_AKTSUB', 'K_P1', 'G12R', 'G12DV', extra=('stage_ERK',))
r.update(pred='P1_primary', cohort='PDAC', note='AKTSUB | ERK matched')
res.append(r)
run('P1', d, ['PDAC'], 'K_P1', 'G12R', 'G12DV', ALLY)
# unadjusted-for-ERK versions of the key PI3K outcomes, plus LSCC has n=1 (skip)
for y in ['stage_MTORC1', 'stage_S6', 'site_PRAS40_T246', 'kin_AKT1']:
    r = C.contrast(d[d.cohort == 'PDAC'], y, 'K_P1', 'G12R', 'G12DV', extra=('stage_ERK',))
    r.update(pred='P1_matched', cohort='PDAC', note='| ERK matched'); res.append(r)
# each vs KRAS-WT is impossible in PDAC (n=4 WT); use G12R vs WT pooled? skip (pre-declared n>=3 rule)

# ---------- P2: KRAS G12C vs G12D / G12V ----------
d2 = D.copy()
for pair, cohs in [(('G12C', 'G12D'), ['LUAD', 'UCEC', 'COAD']),
                   (('G12C', 'G12V'), ['LUAD', 'UCEC']),
                   (('G12V', 'G12D'), ['PDAC', 'COAD', 'UCEC', 'LUAD'])]:
    run(f'P2_{pair[0]}v{pair[1]}', d2, cohs, 'KRAS_allele', pair[0], pair[1], ALLY)
# each allele vs KRAS-WT
for al, cohs in [('G12C', ['LUAD', 'UCEC']), ('G12D', ['COAD', 'UCEC', 'LUAD']),
                 ('G12V', ['COAD', 'UCEC', 'LUAD']), ('G13D', ['UCEC', 'COAD'])]:
    run(f'P2_{al}vWT', d2, cohs, 'KRAS_allele', al, 'WT', ALLY)

# ---------- P3: PIK3CA helical vs kinase ----------
d3 = D.copy()
d3['SGKmAKT'] = d3.site_NDRG1_SGK - d3.site_PRAS40_T246
Y3 = ALLY + ['SGKmAKT', 'kin_SGK3', 'kin_SGK1']
run('P3_helical_v_kinase', d3, ['BRCA', 'UCEC', 'COAD', 'HNSCC'], 'PIK3CA_class', 'helical', 'kinase',
    Y3, extra=('subtype',))
run('P3_helical_v_WT', d3, ['BRCA', 'UCEC', 'COAD', 'LSCC', 'GBM'], 'PIK3CA_class', 'helical', 'WT', Y3, extra=('subtype',))
run('P3_kinase_v_WT', d3, ['BRCA', 'UCEC', 'COAD'], 'PIK3CA_class', 'kinase', 'WT', Y3, extra=('subtype',))
# P3b double vs single hotspot
d3['P3b'] = np.where(d3.PIK3CA_class == 'double', 'double',
                     np.where(d3.PIK3CA_class.isin(['helical', 'kinase']), 'single', 'x'))
run('P3b_double_v_single', d3, ['UCEC', 'BRCA', 'COAD'], 'P3b', 'double', 'single', Y3, extra=('subtype',))
run('P3b_double_v_WT', d3, ['UCEC', 'BRCA'], 'PIK3CA_class', 'double', 'WT', Y3, extra=('subtype',))

# ---------- P4: BRAF V600E vs RAS-mutant / WT in COAD ----------
d4 = D.copy()
rasmut = (d4.KRAS_mut == 1) | (d4.NRAS_mut == 1)
d4['P4'] = np.where(d4.BRAF_class == 'class1', 'BRAFV600',
                    np.where(rasmut & (d4.BRAF_class == 'WT'), 'RASmut',
                             np.where((~rasmut) & (d4.BRAF_class == 'WT'), 'tripleWT', 'x')))
run('P4_V600_v_RASmut', d4, ['COAD'], 'P4', 'BRAFV600', 'RASmut', ALLY, extra=('subtype',))
run('P4_V600_v_tripleWT', d4, ['COAD'], 'P4', 'BRAFV600', 'tripleWT', ALLY, extra=('subtype',))
run('P4_RASmut_v_tripleWT', d4, ['COAD'], 'P4', 'RASmut', 'tripleWT', ALLY, extra=('subtype',))
# class2/3 pooled pan-cancer vs V600 (pre-declared likely underpowered)
d4b = D.copy()
d4b['P4b'] = np.where(d4b.BRAF_class == 'class1', 'class1',
                      np.where(d4b.BRAF_class.isin(['class2', 'class3']), 'class23', 'x'))
d4b['cohort_all'] = 'ALL'
for y in ALLY:
    r = C.contrast(d4b, y, 'P4b', 'class23', 'class1', extra=('cohort',))
    r.update(pred='P4b_class23_v_class1', cohort='PAN', note='cohort-adjusted'); res.append(r)

# ---------- P5: EGFR ex19del vs L858R in LUAD ----------
run('P5_ex19del_v_L858R', D, ['LUAD'], 'EGFR_allele', 'ex19del', 'L858R', ALLY)
run('P5_ex19del_v_WT', D, ['LUAD'], 'EGFR_allele', 'ex19del', 'WT', ALLY)
run('P5_L858R_v_WT', D, ['LUAD'], 'EGFR_allele', 'L858R', 'WT', ALLY)

R = pd.DataFrame(res)
R.to_csv(OUT + 'T1_confirmatory.tsv', sep='\t', index=False)

# sensitivity: same primary contrasts on unadjusted tc
sres = []
d_ = Dtc.copy()
d_['K_P1'] = np.where(d_.KRAS_allele == 'G12R', 'G12R', np.where(d_.KRAS_allele.isin(['G12D', 'G12V']), 'G12DV', 'other'))
d_['SGKmAKT'] = d_.site_NDRG1_SGK - d_.site_PRAS40_T246
rasmut = (d_.KRAS_mut == 1) | (d_.NRAS_mut == 1)
d_['P4'] = np.where(d_.BRAF_class == 'class1', 'BRAFV600', np.where(rasmut & (d_.BRAF_class == 'WT'), 'RASmut', 'x'))
for y in ['stage_AKTSUB', 'stage_MTORC1', 'site_PRAS40_T246']:
    r = C.contrast(d_[d_.cohort == 'PDAC'], y, 'K_P1', 'G12R', 'G12DV', extra=('stage_ERK',)); r.update(pred='P1', cohort='PDAC', variant='tc'); sres.append(r)
for y in ['stage_AKTSUB', 'stage_RTK', 'stage_ADAPT']:
    r = C.contrast(d_[d_.cohort == 'LUAD'], y, 'KRAS_allele', 'G12C', 'G12D'); r.update(pred='P2', cohort='LUAD', variant='tc'); sres.append(r)
    r = C.contrast(d_[d_.cohort == 'COAD'], y, 'P4', 'BRAFV600', 'RASmut'); r.update(pred='P4', cohort='COAD', variant='tc'); sres.append(r)
    r = C.contrast(d_[d_.cohort == 'LUAD'], y, 'EGFR_allele', 'ex19del', 'L858R'); r.update(pred='P5', cohort='LUAD', variant='tc'); sres.append(r)
pd.DataFrame(sres).to_csv(OUT + 'T2_sensitivity_tc.tsv', sep='\t', index=False)

# meta-analyses
metas = []
for (pred, y), g in R.groupby(['pred', 'y']):
    if g.cohort.nunique() < 2: continue
    m = C.meta(g.to_dict('records'))
    if m.get('k', 0) >= 2:
        metas.append({'pred': pred, 'y': y, **m, 'cohorts': ';'.join(g.cohort)})
M = pd.DataFrame(metas)
M.to_csv(OUT + 'T3_meta.tsv', sep='\t', index=False)
print(R.shape, M.shape)
print(R[R.pred.str.startswith('P1')][['pred', 'y', 'n_a', 'n_b', 'beta', 'lo', 'hi', 'p', 'mde']].head(30).to_string())
