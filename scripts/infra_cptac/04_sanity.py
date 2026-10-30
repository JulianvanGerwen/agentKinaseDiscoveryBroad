#!/usr/bin/env python3
"""Sanity checks of CPTAC kinase activity + site-mapping validation. Writes sanity/sanity_checks.tsv and mapping check."""
import os, json
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

OUT = '/workspace/results/infra_cptac'
os.makedirs(f'{OUT}/sanity', exist_ok=True)
COH = ['BRCA', 'CCRCC', 'COAD', 'GBM', 'HNSCC', 'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC']
sam = pd.read_parquet(f'{OUT}/pancan_samples.parquet')
VARS = ['tc', 'tc_padj']
ACT = {(v, m): pd.read_parquet(f'{OUT}/activity/{v}/{m}_z.parquet') for v in VARS + ['paired', 'paired_padj'] for m in ('corev3', 'meanz')}


def cohort_tab(c, name):
    return pd.read_parquet(f'{OUT}/cohorts/{c}/{name}.parquet')


def tumours(c):
    s = sam[(sam.cohort == c) & (sam.sample_type == 'Tumor') & sam.has_phospho]
    return s.set_index('case_id').sample_id


def score(v, m, kin, c):
    A = ACT[(v, m)]
    t = tumours(c)
    if kin not in A.index:
        return pd.Series(np.nan, index=t.index)
    return pd.Series(A.loc[kin].reindex(t.values).values, index=t.index)


rows = []


def test_group(name, kin, cohorts, label_fn, expect='higher'):
    for v in VARS:
        for m in ('corev3', 'meanz'):
            xs, ys = [], []
            for c in cohorts:
                s = score(v, m, kin, c); lab = label_fn(c).reindex(s.index)
                ok = s.notna() & lab.notna()
                xs.append(s[ok]); ys.append(lab[ok].astype(bool))
            x = pd.concat(xs); y = pd.concat(ys)
            if y.sum() < 3 or (~y).sum() < 3:
                rows.append(dict(test=name, kinase=kin, variant=v, method=m, n_pos=int(y.sum()), n_neg=int((~y).sum()), stat=np.nan, p=np.nan, passed=None)); continue
            u, p = mannwhitneyu(x[y], x[~y], alternative='greater' if expect == 'higher' else 'less')
            auc = u / (y.sum() * (~y).sum())
            rows.append(dict(test=name, kinase=kin, variant=v, method=m, n_pos=int(y.sum()), n_neg=int((~y).sum()),
                             stat=round(auc if expect == 'higher' else 1 - auc, 3), p=float(p), passed=bool(p < 0.05)))


def test_corr(name, kin, cohorts, cov_fn):
    for v in VARS:
        for m in ('corev3', 'meanz'):
            xs, ys = [], []
            for c in cohorts:
                s = score(v, m, kin, c); cv = cov_fn(c).reindex(s.index)
                ok = s.notna() & cv.notna()
                xs.append(s[ok]); ys.append(cv[ok] - cv[ok].mean())   # cohort-centred covariate
            x = pd.concat(xs); y = pd.concat(ys)
            if len(x) < 10:
                rows.append(dict(test=name, kinase=kin, variant=v, method=m, n_pos=len(x), stat=np.nan, p=np.nan, passed=None)); continue
            r, p = spearmanr(x, y)
            rows.append(dict(test=name, kinase=kin, variant=v, method=m, n_pos=len(x), n_neg=np.nan, stat=round(r, 3), p=float(p), passed=bool(r > 0 and p < 0.05)))


def gistic(gene, level=2):
    def f(c):
        g = cohort_tab(c, 'cnv_gistic')
        return (g.loc[gene] >= level).astype(float) if gene in g.index else pd.Series(dtype=float)
    return f


def mutated(genes, aa_filter=None):
    def f(c):
        mb = cohort_tab(c, 'mut_binary'); mb.columns = mb.columns.astype(str)
        maf = cohort_tab(c, 'mutations')
        lab = pd.Series(0.0, index=mb.columns)
        for g in genes:
            mm = maf[maf.gene == g]
            if aa_filter:
                mm = mm[mm.protein_change.astype(str).str.contains(aa_filter.get(g, '.'), regex=True)]
            lab[lab.index.isin(mm.case_id)] = 1.0
        return lab
    return f


def mut_or_amp(gene):
    def f(c):
        a = gistic(gene)(c); b = mutated([gene])(c)
        return pd.concat([a, b], axis=1).max(axis=1)
    return f


def annot(col):
    def f(c):
        a = cohort_tab(c, 'annotations').set_index('case_id')
        return a[col] if col in a.columns else pd.Series(dtype=float)
    return f


def rb1_loss(c):
    m = mutated(['RB1'])(c); g = cohort_tab(c, 'cnv_gistic')
    d = (g.loc['RB1'] <= -2).astype(float) if 'RB1' in g.index else pd.Series(dtype=float)
    return pd.concat([m, d], axis=1).max(axis=1)


RAS = {'KRAS': r'p\.G12|p\.G13|p\.Q61', 'NRAS': r'p\.G12|p\.G13|p\.Q61', 'BRAF': r'p\.V600'}
test_group('BRCA ERBB2 amp (GISTIC=2)', 'ERBB2', ['BRCA'], gistic('ERBB2'))
def cbio(col, positive):
    def f(c):
        d = cohort_tab(c, 'subtypes_cbioportal').set_index('case_id')[col].str.lower()
        d = d[d.isin([x.lower() for x in positive]) | d.isin(['negative', 'lumA'.lower(), 'luma', 'lumb', 'basal', 'normal-like'])]
        return d.isin([x.lower() for x in positive]).astype(float)
    return f


test_group('BRCA ERBB2 proteogenomic status Positive (cBioPortal)', 'ERBB2', ['BRCA'], cbio('ERBB2_PROTEOGENOMIC_STATUS', ['Positive']))
test_group('BRCA PAM50 Her2 vs others (cBioPortal)', 'ERBB2', ['BRCA'], cbio('PAM50', ['Her2']))
for k in ('ERBB3', 'EGFR'):   # paralog / shared-substrate check (REPORT.md limitation 3)
    test_group('BRCA ERBB2 amp (GISTIC=2) - ERBB-family paralog', k, ['BRCA'], gistic('ERBB2'))
test_group('GBM EGFR amp or mut', 'EGFR', ['GBM'], mut_or_amp('EGFR'))
test_group('LUAD EGFR mut or amp', 'EGFR', ['LUAD'], mut_or_amp('EGFR'))
for k in ('MAPK1', 'MAPK3'):
    test_group('KRAS/NRAS/BRAF hotspot (LUAD,COAD,UCEC,LSCC,HNSCC; pooled within-cohort)', k, ['LUAD', 'COAD', 'UCEC', 'LSCC', 'HNSCC'], mutated(list(RAS), RAS))
    test_group('KRAS hotspot LUAD', k, ['LUAD'], mutated(['KRAS'], RAS))
for k in ('CDK1', 'CDK2'):
    test_corr('corr with HALLMARK_E2F_TARGETS (all cohorts, cohort-centred)', k, COH, annot('HALLMARK_E2F_TARGETS'))
    test_corr('corr with HALLMARK_G2M_CHECKPOINT (all cohorts)', k, COH, annot('HALLMARK_G2M_CHECKPOINT'))
test_group('RB1 mut or deep del (GBM,LSCC,BRCA,UCEC pooled)', 'CDK2', ['GBM', 'LSCC', 'BRCA', 'UCEC', 'OV'], rb1_loss)
test_group('PIK3CA hotspot/PTEN mut (BRCA,UCEC)', 'AKT1', ['BRCA', 'UCEC'], mutated(['PIK3CA', 'PTEN'], {'PIK3CA': r'p\.E542|p\.E545|p\.H1047', 'PTEN': '.'}))
test_corr('corr with BCM PTM-SEA ERK2/MAPK1 score (same data; agreement)', 'MAPK1', COH, annot('PTM_SEA_KINASE-PSP_ERK2/MAPK1'))
test_corr('corr with BCM PTM-SEA CDK1 score (agreement)', 'CDK1', COH, annot('PTM_SEA_KINASE-PSP_CDK1'))
test_corr('corr with BCM PTM-SEA Akt1 score (agreement)', 'AKT1', COH, annot('PTM_SEA_KINASE-PSP_Akt1/AKT1'))
test_corr('corr with PROGENy_MAPK (RNA)', 'MAPK1', COH, annot('PROGENy_MAPK'))
test_corr('corr with PROGENy_EGFR (RNA)', 'EGFR', COH, annot('PROGENy_EGFR'))

# tumour vs paired normal: CDK1/CDK2 up in tumours (one-sample sign test on paired z)
from scipy.stats import wilcoxon
for v in ('paired', 'paired_padj'):
    for m in ('corev3', 'meanz'):
        A = ACT[(v, m)]
        for k in ('CDK1', 'CDK2', 'AURKA', 'PLK1'):
            if k not in A.index: continue
            x = A.loc[k].dropna()
            if len(x) < 10: continue
            st, p = wilcoxon(x, alternative='greater')
            rows.append(dict(test='tumour > paired normal (all paired cases)', kinase=k, variant=v, method=m, n_pos=len(x),
                             stat=round(float((x > 0).mean()), 3), p=float(p), passed=bool(p < 0.05)))

res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/sanity/sanity_checks.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
print(res.to_string())

# ---- corev3 vs meanz agreement (per kinase across tumours, tc)
a, b = ACT[('tc', 'corev3')], ACT[('tc', 'meanz')]
T = sam[(sam.sample_type == 'Tumor') & sam.has_phospho].sample_id
ks = [k for k in b.index if k in a.index]
rs = {}
for k in ks:
    x, y = a.loc[k, T], b.loc[k, T]; ok = x.notna() & y.notna()
    if ok.sum() > 50: rs[k] = spearmanr(x[ok], y[ok])[0]
rs = pd.Series(rs)
agree = {'n_kinases': int(len(rs)), 'median_spearman': float(rs.median()), 'frac_gt_0.5': float((rs > 0.5).mean())}
# ---- mapping validation vs PSP human site list
ps = pd.read_csv('/workspace/data/Phosphorylation_site_dataset', sep='\t', skiprows=3, dtype=str, encoding='latin-1')
ps = ps[ps.ORGANISM.eq('human')]
pk = set(ps.ACC_ID + '_' + ps.MOD_RSD.str.replace('-p', '', regex=False))
smap = pd.read_parquet(f'{OUT}/mapping/site_map.parquet')
mk = smap.site_key.dropna().unique()
# negative control: keys with GENCODE position kept but residue shifted by +1 should rarely be in PSP
shifted = smap.dropna(subset=['site_key'])
neg = (shifted.uniprot + '_' + shifted.residue + (shifted.pos_uniprot + 2).astype(str)).unique()
ch = smap[smap.method != 'same_pos'].site_key.dropna().unique()
mv = {'mapped_keys': len(mk), 'frac_in_PSP_human': float(np.mean([k in pk for k in mk])),
      'frac_in_PSP_for_repositioned_sites': float(np.mean([k in pk for k in ch])), 'n_repositioned': len(ch),
      'neg_control_pos+2_frac_in_PSP': float(np.mean([k in pk for k in neg]))}
json.dump({'corev3_vs_meanz_tc': agree, 'mapping_vs_PSP': mv}, open(f'{OUT}/sanity/agreement_and_mapping.json', 'w'), indent=1)
print(agree, mv)
