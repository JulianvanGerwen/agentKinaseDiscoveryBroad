# infra_cptac — CPTAC pan-cancer kinase-signalling assets (built 2026-10-29)

**What**: 10 CPTAC cohorts (BCM GENCODE-v34 harmonised freeze v1.2, Li et al. 2023), with phosphosites keyed to canonical UniProt, plus protein, RNA, mutations, CNV, clinical, immune and subtype tables and per-sample kinase activity. Scripts are in `/workspace/scripts/infra_cptac/` (00 download, 01 mapping, 02 cohorts, 03 activity, 04 sanity, 05 subtypes, `cptac_load.py`). Python env: `/home/plm_env/bin/python`. The whole pipeline (02→04) re-runs in about 2 min.

**Sources**: phospho, protein and BCM meta come from the local read-only `kinaseSubMain/external/cptac_pancan_v1.2` and `cptac_pancancer_pdc`. I downloaded the following myself (they were missing locally): PDC `Mutation_BCM_v1`, `CNV_BCM_v1` and `RNA_BCM_v1` (into `/workspace/cache/cptac/`), and cohort-paper clinical/subtype attributes from cBioPortal.

## Layout (`/workspace/results/infra_cptac/`)
| path | schema |
|---|---|
| `cohorts/<C>/phospho.parquet` | site_key (`P04626_Y1248`) × sample, float32 log2 ref-intensity; tumour + normal columns |
| `cohorts/<C>/protein.parquet` | UniProt × sample (float32). `rna.parquet`: gene symbol × sample (RSEM UQ log2) |
| `cohorts/<C>/cnv_log2.parquet`, `cnv_gistic.parquet` | gene × case (WES log2 ratio; GISTIC −2..2) |
| `cohorts/<C>/mut_binary.parquet`, `mutations.parquet` | gene × case (BCM non-synonymous); long MAF with gene, protein_change, class, VAF |
| `cohorts/<C>/clinical.parquet`, `annotations.parquet` | per case: age, sex, stage, grade, OS/PFS days+event, purity, TMB, MSI, POLE, CIN; CIBERSORT/xCell/ESTIMATE/PROGENy/Hallmark ssGSEA/PTM-SEA/mutational signatures |
| `cohorts/<C>/subtypes_cbioportal.parquet` | cBioPortal paper studies: BRCA PAM50, NMF and HER2/ER status; LUAD/LSCC NMF and immune clusters; UCEC/COAD MSI and genomic subtype; GBM and PDAC clinical (strings) |
| `cohorts/<C>/input_{tc,paired,tc_padj,paired_padj}.parquet` | the site × sample statistics used for activity |
| `pancan_samples.parquet` | sample_id (`<case>_T` / `<case>_N`), case_id, cohort, sample_type, has_{phospho,protein,rna,mut,cnv,clinical,paired_normal_phospho} |
| `pancan_clinical.parquet`, `cohort_counts.tsv` | stacked clinical; per-cohort counts |
| `activity/<variant>/{corev3_z,corev3_zcore,corev3_n,meanz_z,meanz_n}.parquet` | kinase × sample, all cohorts (1,553 samples for tc/tc_padj; 492 for the paired variants) |
| `mapping/site_map.parquet` | CPTAC ID → uniprot, pos, method, residue_check |
| `edges/`, `model_corev3/` | the edge library and core-v3 config I used |
| `sanity/sanity_checks.tsv` | sanity-check results |

**Counts (tumour/normal with phospho)**: BRCA 122/0 · CCRCC 103/80 · COAD 97/100 · GBM 99/0 · HNSCC 108/62 · LSCC 108/99 · LUAD 110/101 · OV 83/19 · PDAC 105/44 · UCEC 95/18. That is 1,030 tumours and 523 normals. Paired T–N cases: 492; the UCEC normals are unpaired NX donors. Mapped sites per cohort range from 35k (COAD) to 66k (UCEC). Mutation, CNV and RNA data exist for about 100% of tumours; OS is available for 944 tumours.

## Site mapping
- **Mapping rate: 99.70%** (126,659 of 127,040 IDs, giving 125,753 unique keys).
- **Method**: ENSP → SwissProt (BCM mapping table), with gene-name fallback for 4.3% of IDs. The site's 15-mer is then located in the canonical sequence. 97.2% of sites keep the same position; 2.7% are repositioned by window match.
- **Residue check**: the residue letter matches the canonical sequence for 100% of mapped sites.
- **Agreement with PSP**: 61% of mapped keys are human PSP sites. A negative control with the position shifted by +2 gives 6.6%.

## Activity (inputs and methods)
**Input variants**
- `tc`: per-site log2 minus the cohort's tumour median, applied to tumours and normals. Only sites seen in at least 25% of tumours and at least 10 tumours are used (17–41k sites per cohort).
- `paired`: tumour − matched normal.
- `tc_padj`: the residual of phospho ~ a + b·protein (same UniProt). The fit is OLS per site, on tumours, needing at least 10 pairs; median slope b is 0.43–0.84. Sites without protein (about 7%) are dropped.
- `paired_padj`: the padj residual for the tumour minus the padj residual for the normal.

**Methods**
- `corev3`: the kinaseActivityMain MVP core-v3, scoring 503 kinases. It reports z, z_core and n_substrates.
- `meanz`: the per-sample mean robust z of PSP-curated substrates, requiring at least 5 measured substrates (about 100 kinases per sample).
- The two methods agree per kinase across tumours with a median Spearman of 0.53 (130 kinases).

## Sanity checks
AUROC or Spearman ρ; values are for `tc`, with `tc_padj` in brackets where noted.

| check | corev3 | meanz | verdict |
|---|---|---|---|
| MAPK1/3 vs KRAS/NRAS/BRAF hotspot (5 cohorts, n=124/390) | 0.58 / 0.60 | 0.56 / 0.56 | PASS (all variants) |
| MAPK1/3 vs KRAS hotspot, LUAD | 0.62 / 0.65 | 0.62 / 0.60 | PASS |
| EGFR vs amp/mut, GBM (48/51) | 0.83 [padj 0.55] | 0.84 [0.60] | PASS (tc); padj weakens it, as expected for autophosphorylation sites |
| EGFR vs mut/amp, LUAD | 0.69 [0.62] | 0.57 | PASS corev3 / FAIL meanz |
| AKT1 vs PIK3CA hotspot/PTEN mut (BRCA+UCEC) | 0.62 | 0.65 | PASS |
| CDK1/CDK2 vs E2F/G2M Hallmark (n=1027) | ρ 0.63–0.73 [0.48–0.54] | ρ 0.68–0.76 | PASS |
| CDK2 vs RB1 mut/deep deletion (n=33) | 0.63 [0.58 ns] | 0.61 | PASS (tc) |
| CDK1/CDK2 tumour > paired normal (n=492) | 86% / 95% of pairs > 0 | 90% / 91% | PASS |
| **ERBB2 vs ERBB2 amp, BRCA (16/97)** | **0.50** (PAM50-Her2: 0.68, p=0.02) | not scorable | **FAIL**: ERBB3 picks the amp up instead (0.79, p=1e-4) |
| EGFR vs PROGENy-EGFR (RNA) | ρ ≈ 0 | ρ ≈ 0 | FAIL |

## Caveats
- **ERBB2**: only about 7 measured substrates, and its signal is assigned to ERBB3/EGFR. This is the paralog problem noted in REPORT.md. Use kinase families for paralogs, and use ERBB2 protein or pY1248 directly.
- **Calibration**: |z|>1.96 occurs for 15–19% of kinase-sample pairs, because tumours really are heterogeneous. Treat z as a ranking, not a p-value.
- **Protein adjustment**:
  - It removes the protein-abundance component, including real autophosphorylation and amplification signal (see GBM EGFR).
  - Normals are adjusted with the tumour fit.
  - The tumour median for `tc` mixes subtypes, so z is relative to the cohort.
- **Missing values**: −Inf in the source files (log2 of 0) is set to NaN. Nothing is imputed.
- **Sample coverage**: BRCA and GBM have no normals. OV has only 10 paired cases.
- **Mutation coordinates**: MAF protein changes use BCM transcript coordinates, not UniProt.
- **Subtypes**: none from cBioPortal for CCRCC, HNSCC or OV. These are in the cohort-paper supplements (Clark 2019, Huang 2021, McDermott 2020) or on LinkedOmics.
- **Edge library**: mine is `edges/edges_psp_omnipath.tsv.gz`, built with `mvp/build_edges.py` from /workspace/data PSP + OmniPath; it is identical to `mvp/model/edges.tsv.gz`. infra_compendium uses the research `curated_union` instead; the two share 25,128 curated edges.

## Loading
```python
import sys; sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
ph = cl.matrix('LUAD', 'phospho'); sam = cl.samples(); z = cl.activity('tc', 'corev3_z')[cl.tumour_ids('LUAD')]
mut = cl.table('LUAD', 'mutations'); clin = cl.table('LUAD', 'clinical'); n = cl.activity('tc_padj', 'corev3_n')
```
