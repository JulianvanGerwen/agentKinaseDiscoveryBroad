#!/usr/bin/env python3
"""Build analysis-ready per-cohort CPTAC assets (BCM-harmonised freeze v1.2 + PDC pan-cancer files).

Outputs under /workspace/results/infra_cptac/cohorts/<COHORT>/ :
  phospho.parquet       site_key x sample (float32, log2 ref-intensity; tumour + normal samples)
  protein.parquet       UniProt x sample (float32, log2)
  rna.parquet           gene symbol x sample (float32, log2 RSEM UQ; tumour [+ normal])
  cnv_log2.parquet      gene symbol x case (float32, WES log2 ratio);  cnv_gistic.parquet (int8, -2..2)
  mut_binary.parquet    gene symbol x case (int8; BCM non-synonymous gene-level calls)
  mutations.parquet     long MAF subset (case, gene, variant class, protein change, VAF), non-silent coding only
  samples.parquet       sample_id, case_id, cohort, sample_type (Tumor/Normal), has_* flags
  clinical.parquet      per case: age, sex, stage, grade, OS/PFS, purity, TMB, MSI, ...
  annotations.parquet   per case: immune deconvolution / ESTIMATE / PROGENy / Hallmark / PTM-SEA / mut. signatures
Sample IDs: '<case_id>_T' (tumour) and '<case_id>_N' (normal adjacent tissue).
"""
import os, json, re
import numpy as np, pandas as pd

SRC = '<HACKATHON_ROOT>/kinaseSubMain/external'
V12 = f'{SRC}/cptac_pancan_v1.2'
CACHE = '/workspace/cache/cptac'
RAW = f'{CACHE}/raw'
OUT = '/workspace/results/infra_cptac'
COH = ['BRCA', 'CCRCC', 'COAD', 'GBM', 'HNSCC', 'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC']
ALT = {'OV': 'HGSC', 'CCRCC': 'ccRCC'}            # naming in PDC RNA/mutation/CNV/clinical files
MAPTAB = f'{CACHE}/bcm_readme_phospho/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt'

smap = pd.read_parquet(f'{OUT}/mapping/site_map.parquet').set_index('cptac_id')
mt = pd.read_csv(MAPTAB, sep='\t', usecols=['gene', 'gene_name', 'SwissProt', 'Primary_select'], dtype=str)
g2name = mt.drop_duplicates('gene').set_index('gene')['gene_name']
g2sp = mt[mt.SwissProt.notna()].sort_values('Primary_select', ascending=False).drop_duplicates('gene').set_index('gene')['SwissProt']


def read_tab(path):
    df = pd.read_csv(path, sep='\t', index_col=0, low_memory=False)
    df = df.replace([np.inf, -np.inf], np.nan)
    df.index.name = None
    return df


def collapse(df, key):
    """Re-key rows by `key` (Series aligned to df.index); drop NaN keys; duplicates -> row with most observations."""
    k = key.reindex(df.index)
    df = df[k.notna().values]
    k = k[k.notna()]
    nobs = df.notna().sum(axis=1).values
    order = np.lexsort((-nobs, k.values))
    df = df.iloc[order]; k = k.iloc[order]
    keep = ~pd.Index(k.values).duplicated()
    df = df[keep]; df.index = k.values[keep]
    return df.sort_index()


def f32(df):
    return df.astype(np.float32)


def omic_pair(cohort, kind):
    """Tumour + normal matrices for phospho/proteome from the BCM v1.2 freeze."""
    stem = {'phospho': 'phospho_site_abundance', 'protein': 'proteomics_gene_abundance'}[kind]
    out = []
    for st, suf in (('Tumor', '_T'), ('Normal', '_N')):
        p = f'{V12}/{cohort}/{cohort}_{stem}_log2_reference_intensity_normalized_{st}.txt'
        if os.path.exists(p):
            d = read_tab(p); assert d.columns.is_unique
            d = d.replace([np.inf, -np.inf], np.nan)          # source files contain -Inf (log2 of 0 intensity) -> missing
            d.columns = [c + suf for c in d.columns]; out.append(d)
    return pd.concat(out, axis=1)


def main():
    stats = {}
    allsamples, allclin = [], []
    for c in COH:
        a = ALT.get(c, c)
        od = f'{OUT}/cohorts/{c}'; os.makedirs(od, exist_ok=True)
        s = {}
        # ---- phospho
        ph = omic_pair(c, 'phospho')
        s['phospho_ids'] = len(ph)
        ph = collapse(ph, smap['site_key'])
        s['phospho_sites_mapped'] = len(ph)
        f32(ph).to_parquet(f'{od}/phospho.parquet')
        # ---- protein (ENSG -> SwissProt canonical)
        pr = omic_pair(c, 'protein')
        s['protein_genes'] = len(pr)
        prot_ann = pd.DataFrame({'ensg': pr.index, 'gene': g2name.reindex(pr.index).values, 'uniprot': g2sp.reindex(pr.index).values})
        pr = collapse(pr, g2sp)
        s['protein_uniprot'] = len(pr)
        f32(pr).to_parquet(f'{od}/protein.parquet')
        # ---- RNA
        rna = []
        for st, suf in (('Tumor', '_T'), ('Normal', '_N')):
            p = f'{RAW}/RNA_BCM_v1/{a}_RNAseq_gene_RSEM_coding_UQ_1500_log2_{st}.txt'
            if os.path.exists(p):
                d = read_tab(p); d.columns = [x + suf for x in d.columns]; rna.append(d)
        if rna:
            rna = pd.concat(rna, axis=1); rna = collapse(rna, g2name)
            f32(rna).to_parquet(f'{od}/rna.parquet'); s['rna_genes'] = len(rna); s['rna_samples'] = rna.shape[1]
        else:
            rna = None
        # ---- CNV
        cl = collapse(read_tab(f'{RAW}/CNV_BCM_v1/{a}_WES_CNV_gene_ratio_log2.txt'), g2name)
        f32(cl).to_parquet(f'{od}/cnv_log2.parquet')
        cg = collapse(read_tab(f'{RAW}/CNV_BCM_v1/{a}_WES_CNV_gene_gistic_level.txt'), g2name)
        cg.astype('float32').round().astype('Int8').to_parquet(f'{od}/cnv_gistic.parquet')
        s['cnv_cases'] = cl.shape[1]
        # ---- mutations
        mb = collapse(read_tab(f'{RAW}/Mutation_BCM_v1/{a}_somatic_mutation_gene_level_binary.txt'), g2name)
        mb.astype('int8').to_parquet(f'{od}/mut_binary.parquet'); s['mut_cases'] = mb.shape[1]
        maf = pd.read_csv(f'{RAW}/Mutation_BCM_v1/{a}_somatic_mutation.maf', sep='\t', low_memory=False,
                          usecols=['Hugo_Symbol', 'gene_name', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                                   'Variant_Classification', 'Variant_Type', 'tx', 'aaChange', 'Tumor_Sample_Barcode', 'TumorVAF', 't_alt_count', 't_ref_count'])
        nonsil = {'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',
                  'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site'}
        maf = maf[maf.Variant_Classification.isin(nonsil)].rename(columns={'Tumor_Sample_Barcode': 'case_id', 'Hugo_Symbol': 'gene',
                                                                         'gene_name': 'ensg', 'aaChange': 'protein_change'})
        maf['case_id'] = maf.case_id.astype(str)
        maf.reset_index(drop=True).to_parquet(f'{od}/mutations.parquet'); s['maf_nonsilent'] = len(maf)
        # ---- clinical / annotations
        meta = pd.read_csv(f'{CACHE}/clinical/Clinical_meta_data_v1/{a}_meta.txt', sep='\t', index_col=0, low_memory=False, dtype=str)
        meta = meta.drop(index='data_type'); meta.index = meta.index.astype(str); meta.index.name = 'case_id'
        ann_pat = re.compile(r'^(CIBERSORT_|xCell_|ESTIMATE_|PROGENy_|HALLMARK_|PTM_SEA_|Mutation_signature_|CNV_index_chr)')
        acols = [x for x in meta.columns if ann_pat.match(x)]
        ccols = [x for x in meta.columns if x not in acols]
        clin = meta[ccols].copy()
        for col in clin.columns:
            num = pd.to_numeric(clin[col], errors='coerce')
            if num.notna().sum() >= 0.9 * clin[col].notna().sum() and clin[col].notna().any() and col not in ('Tumor', 'Normal'):
                clin[col] = num
        clin.insert(0, 'cohort', c)
        clin.reset_index().to_parquet(f'{od}/clinical.parquet', index=False)
        ann = meta[acols].apply(pd.to_numeric, errors='coerce').astype('float32'); ann.insert(0, 'cohort', c)
        ann.reset_index().to_parquet(f'{od}/annotations.parquet', index=False)
        # ---- sample table
        sam = pd.DataFrame({'sample_id': sorted(set(ph.columns) | set(pr.columns))})
        sam['case_id'] = sam.sample_id.str[:-2]; sam['sample_type'] = np.where(sam.sample_id.str.endswith('_T'), 'Tumor', 'Normal')
        sam.insert(2, 'cohort', c)
        sam['has_phospho'] = sam.sample_id.isin(ph.columns); sam['has_protein'] = sam.sample_id.isin(pr.columns)
        sam['has_rna'] = sam.sample_id.isin(rna.columns) if rna is not None else False
        sam['has_mut'] = sam.case_id.isin(mb.columns.astype(str)); sam['has_cnv'] = sam.case_id.isin(cl.columns.astype(str))
        sam['has_clinical'] = sam.case_id.isin(clin.index)
        sam['has_paired_normal_phospho'] = sam.case_id.isin(sam.case_id[(sam.sample_type == 'Normal') & sam.has_phospho]) & \
            sam.case_id.isin(sam.case_id[(sam.sample_type == 'Tumor') & sam.has_phospho])
        sam.to_parquet(f'{od}/samples.parquet', index=False)
        allsamples.append(sam)
        allclin.append(clin.reset_index())
        s.update({'tumor_phospho': int(((sam.sample_type == 'Tumor') & sam.has_phospho).sum()),
                  'normal_phospho': int(((sam.sample_type == 'Normal') & sam.has_phospho).sum()),
                  'paired_TN_phospho_cases': int(sam[sam.has_paired_normal_phospho].case_id.nunique()),
                  'tumor_with_mut': int(((sam.sample_type == 'Tumor') & sam.has_mut & sam.has_phospho).sum()),
                  'tumor_with_rna': int(((sam.sample_type == 'Tumor') & sam.has_rna & sam.has_phospho).sum()),
                  'tumor_with_OS': int(clin.reindex(sam[sam.sample_type == 'Tumor'].case_id).OS_days.notna().sum())})
        stats[c] = s
        print(c, s, flush=True)
    allsam = pd.concat(allsamples, ignore_index=True)
    allsam.to_parquet(f'{OUT}/pancan_samples.parquet', index=False)
    pd.concat(allclin, ignore_index=True).to_parquet(f'{OUT}/pancan_clinical.parquet', index=False)
    pd.DataFrame(stats).T.to_csv(f'{OUT}/cohort_counts.tsv', sep='\t')


if __name__ == '__main__':
    main()
