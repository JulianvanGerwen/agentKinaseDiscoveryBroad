"""I29 step 1: assemble per-tumour LUAD table (groups, covariates, ischaemia, stage outcomes, Broad ERK peptides, RNA, protein).
Output: /workspace/cache/b_I29/luad.parquet (+ LSCC etc. in 04)."""
import sys, gzip, re, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
sys.path.insert(0, '/workspace/scripts/audit_I10'); import common as c10
CA = '/workspace/cache/b_I29/'
PDC = '<HACKATHON_ROOT>/kinaseSubMain/external/cptac_pancancer_pdc/'

def zs(s): return (s - s.mean()) / s.std()

# ---------------- groups / covariates (all cohorts used later) ----------------
AL = pd.read_parquet('/workspace/results/b_I27/tables/alleles.parquet')

def build(coh):
    a = AL[AL.cohort == coh].set_index('case_id').copy()
    ras = a[['KRAS_mut', 'NRAS_mut', 'BRAF_mut', 'NF1_mut']].sum(axis=1) > 0
    egfr = a.EGFR_mut == 1
    a['group'] = np.where(egfr, 'EGFR', np.where(ras, 'RAS', 'WT'))
    try:
        sb = cl.table(coh, 'subtypes_cbioportal').drop_duplicates('case_id').set_index('case_id')
    except Exception:
        sb = pd.DataFrame(index=a.index)
    gis = cl.matrix(coh, 'cnv_gistic')
    a['EGFR_gistic'] = gis.loc['EGFR'].reindex(a.index).values if 'EGFR' in gis.index else np.nan
    for c in ['COUNTRY_OF_ORIGIN', 'SMOKING_STATUS', 'PUTATIVE_DRIVER_MUTATION', 'ALK_FUSION', 'RET_FUSION', 'ROS1_FUSION', 'EXPERIMENT', 'STAGE', 'SEX']:
        a[c] = sb[c].reindex(a.index).values if c in sb.columns else np.nan
    clin = cl.table(coh, 'clinical').drop_duplicates('case_id').set_index('case_id')
    for c in ['Sex', 'Stage', 'Tobacco_smoking_history', 'Age']:
        a['clin_' + c] = clin[c].reindex(a.index).values if c in clin.columns else np.nan
    fus = (a[['ALK_FUSION', 'RET_FUSION', 'ROS1_FUSION']].astype(str) == '1').any(axis=1)
    other_drv = a.PUTATIVE_DRIVER_MUTATION.astype(str).str.contains('EGFR-AMP|ERBB2|HRAS|NRAS|BRAF', regex=True)
    a['WT_clean'] = (a.group == 'WT') & ~fus & ~other_drv & ~(a.EGFR_gistic == 2)
    a['EGFR_clean'] = (a.group == 'EGFR') & a.EGFR_allele.isin(['ex19del', 'L858R', 'multi']) & ~ras
    a['RAS_any'] = ras
    return a

# ---------------- outcomes ----------------
def outcomes(variant):
    O = pd.read_parquet(f'/workspace/results/b_I27/tables/outcomes_{variant}.parquet')
    keep = [c for c in O.columns if c.startswith(('stage_', 'fp_', 'site_'))] + ['purity', 'prolif', 'cohort', 'case_id']
    return O[keep].set_index('case_id')

def isch_scores(coh):
    """own ISCH and ISCH_noERK (tc_padj), plus val_E isch/isch_noPath."""
    ISCH = c10.ischaemia_set(200)
    psp = pd.read_csv('/workspace/data/Kinase_Substrate_Dataset', sep='\t', skiprows=3, encoding='latin-1')
    psp = psp[(psp.KIN_ORGANISM == 'human') & (psp.SUB_ORGANISM == 'human')]
    erksub = set((psp[psp.KIN_ACC_ID.isin(['P28482', 'P27361'])].SUB_ACC_ID + '_' + psp[psp.KIN_ACC_ID.isin(['P28482', 'P27361'])].SUB_MOD_RSD))
    FP = c10.FP
    fpsites = set(FP[FP.group.isin(['KI_MEK', 'KI_ERK', 'KI_RAF', 'EGF', 'GF_all', 'MAPKaxis'])].site)
    prot_excl = {'P27361', 'P28482', 'Q15418', 'P51812'}
    keep = [s for s in ISCH.index if s not in erksub and s not in fpsites and s.split('_')[0] not in prot_excl]
    Z = c10.standardise(c10.cohort(coh, 'input_tc_padj'))
    d = pd.DataFrame({'ISCH': c10.score(Z, ISCH), 'ISCH_noERK': c10.score(Z, ISCH[keep])})
    d['n_isch_noERK_sites_def'] = len(keep)
    d.index = [i[:-2] for i in d.index]
    v = pd.read_csv('/workspace/results/val_E/ischaemia_scores.tsv', sep='\t', index_col=0)
    v = v[(v.cohort == coh) & (v.variant == 'input_tc_padj')].set_index('case_id')
    d['isch_valE'] = v.isch.reindex(d.index); d['isch_noPath_valE'] = v.isch_noPath.reindex(d.index)
    return d, len(ISCH), len(keep)

# ---------------- Broad peptide-level ERK ----------------
PEP = {'MAPK1': 'NP_002736.3', 'MAPK3': 'NP_002737.2'}
PROT_RS = {'MAPK1': 'NP_002736.3', 'MAPK3': 'NP_002737.2', 'DUSP4': 'NP_001385.1', 'DUSP6': 'NP_001937.2'}
def norm_col(c):
    c = c.replace('.', '-')
    c = re.sub(r'-T$', '', c)
    if re.match(r'^X\d', c): c = c[1:]
    return c

def broad(cases):
    rows = {}
    with gzip.open(PDC + 'Phosphoproteome_Broad_Institute_harmonized_v1.tsv.gz', 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            if line.startswith(('NP_002736.3_', 'NP_002737.2_')):
                p = line.rstrip('\n').split('\t'); rows[p[0]] = p[1:]
    cols = [norm_col(c) for c in hdr[1:]]
    P = pd.DataFrame(rows, index=cols).replace('', np.nan).astype(float)
    prow = {}
    with gzip.open(PDC + 'Proteome_Broad_Institute_harmonized_v1.tsv.gz', 'rt') as f:
        hdr2 = f.readline().rstrip('\n').split('\t')
        want = set(PROT_RS.values())
        for line in f:
            k = line.split('\t', 1)[0]
            if k in want:
                p = line.rstrip('\n').split('\t'); prow[k] = p[1:]
    Q = pd.DataFrame(prow, index=[norm_col(c) for c in hdr2[1:]]).replace('', np.nan).astype(float)
    Q.columns = [{v: k for k, v in PROT_RS.items()}[c] + '_protBroad' for c in Q.columns]
    P = P[~P.index.duplicated()].reindex(cases); Q = Q[~Q.index.duplicated()].reindex(cases)
    out = pd.DataFrame(index=cases)
    for g, rs in PEP.items():
        for tag, pat in [('pT', r'_T(185|202)t_1_1'), ('pY', r'_Y(187|204)y_1_1'), ('pTpY', r'_T(185|202)tY(187|204)y_2_2')]:
            col = [c for c in P.columns if c.startswith(rs) and re.search(pat, c)]
            out[f'{g}_{tag}_raw'] = P[col[0]] if col else np.nan
            # protein-adjusted: residual on Broad protein of same gene (OLS on available pairs)
            y = out[f'{g}_{tag}_raw']; x = Q[g + '_protBroad']; m = y.notna() & x.notna()
            if m.sum() >= 10:
                b = np.polyfit(x[m], y[m], 1); out[f'{g}_{tag}'] = y - (b[0] * x + b[1])
            else:
                out[f'{g}_{tag}'] = np.nan
    return out.join(Q)

GENES = ['ETV4', 'ETV5', 'ETV1', 'SPRED1', 'SPRED2', 'DUSP6', 'DUSP4', 'DUSP5', 'SPRY2', 'SPRY4', 'PHLDA1', 'EGR1', 'FOS', 'FOSB', 'JUN', 'IER3', 'SOCS3', 'MAPK1', 'MAPK3', 'EGFR']
PROT = {'MAPK1': 'P28482', 'MAPK3': 'P27361', 'DUSP6': 'Q16828', 'DUSP4': 'Q13115', 'SPRY2': 'O43597', 'SPRY4': 'Q9C004',
        'SPRED1': 'Q7Z699', 'SPRED2': 'Q7Z698', 'ETV4': 'P43268', 'ETV5': 'P41161', 'EGFR': 'P00533', 'DUSP5': 'Q16690'}

def assemble(coh, with_broad=True):
    a = build(coh)
    tum = cl.tumour_ids(coh); cases = [t[:-2] for t in tum]
    a = a.reindex(cases)
    for v in ['tc_padj', 'tc']:
        O = outcomes(v); O = O[O.cohort == coh].drop(columns='cohort')
        suf = '' if v == 'tc_padj' else '_tc'
        O = O.add_suffix(suf).rename(columns={'purity' + suf: 'purity', 'prolif' + suf: 'prolif'})
        if v == 'tc': O = O.drop(columns=['purity', 'prolif'])
        a = a.join(O)
    rna = cl.matrix(coh, 'rna'); rna.columns = [c[:-2] if c.endswith('_T') else c + '__N' for c in rna.columns]
    for g in GENES:
        a['rna_' + g] = rna.loc[g].reindex(a.index).values if g in rna.index else np.nan
    pr = cl.matrix(coh, 'protein'); pr = pr[[c for c in pr.columns if c.endswith('_T')]]; pr.columns = [c[:-2] for c in pr.columns]
    for g, u in PROT.items():
        a['prot_' + g] = pr.loc[u].reindex(a.index).values if u in pr.index else np.nan
    ann = cl.table(coh, 'annotations').drop_duplicates('case_id').set_index('case_id')
    a['PROGENy_MAPK'] = ann['PROGENy_MAPK'].reindex(a.index).values if 'PROGENy_MAPK' in ann else np.nan
    a['PROGENy_EGFR'] = ann['PROGENy_EGFR'].reindex(a.index).values if 'PROGENy_EGFR' in ann else np.nan
    d, n1, n2 = isch_scores(coh); a = a.join(d)
    t = pd.read_csv(CA + 'gdc_ischaemia_times.tsv', sep='\t').set_index('case_id')
    a['excision_min'] = t.excision_min.reindex(a.index); a['clamping_min'] = t.clamping_min.reindex(a.index)
    if with_broad: a = a.join(broad(list(a.index)))
    a['cohort'] = coh
    print(coh, a.shape, 'ISCH sites', n1, 'noERK', n2, a.group.value_counts().to_dict(), flush=True)
    return a

if __name__ == '__main__':
    for coh in ['LUAD', 'LSCC', 'HNSCC', 'GBM']:
        A = assemble(coh); A.to_parquet(CA + f'{coh}.parquet')
        if coh in ('LUAD', 'LSCC'):
            b = [c for c in A.columns if c.startswith('MAPK') and ('_p' in c)]
            print(A[b].notna().sum().to_dict())
