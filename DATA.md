
## CPTAC pan-cancer (infra_cptac, 2026-10-29)
- **Assets**: `/workspace/results/infra_cptac/`. Read `SUMMARY.md` there first. Loader: `/workspace/scripts/infra_cptac/cptac_load.py`. Env: `/home/plm_env/bin/python`.
- **Cohorts**: 10 (BRCA CCRCC COAD GBM HNSCC LSCC LUAD OV PDAC UCEC), with 1,030 tumours and 523 normal adjacent tissue samples with phospho; 492 tumour–normal pairs. Source: BCM GENCODE-v34 harmonised freeze v1.2 (Li et al. 2023 *Cancer Cell*), read-only at `kinaseSubMain/external/cptac_pancan_v1.2` and `cptac_pancancer_pdc`.
- **Per cohort** (`cohorts/<C>/`):
  - phospho (UniProt site_key × sample; 99.7% of sites mapped, residue-validated) and protein (UniProt × sample), both float32 parquet
  - RNA, CNV (log2/GISTIC), mutations (binary + long MAF with protein change)
  - clinical (OS/PFS, stage, age, sex, purity, TMB), annotations (immune deconvolution, PROGENy, Hallmark, PTM-SEA) and cBioPortal subtypes (BRCA PAM50, LUAD/LSCC NMF, UCEC MSI/genomic)
- **Kinase activity** (`activity/<variant>/`; kinase × sample):
  - variants: `tc` (centred on the cohort tumour median), `paired` (T−N), `tc_padj` and `paired_padj` (phospho residualised on protein)
  - methods: core-v3 (z, n_substrates) and meanz (mean z of PSP substrates, at least 5 substrates)
- **Downloaded by me**: PDC Mutation/CNV/RNA_BCM_v1 zips (`/workspace/cache/cptac/pdc_dl`) and cBioPortal clinical JSON (`/workspace/cache/cptac/cbioportal`).
- **Not available**: subtypes for CCRCC, HNSCC and OV (they are in the cohort-paper supplements). Proteome/phospho from the UMich or Broad pipelines are on disk but were not processed.

## Phosphoproteomics compendium (infra_compendium, 2026-10-29)
- **Assets**: `/workspace/results/infra_compendium/`. Read `SUMMARY.md` there first. Loader: `/workspace/scripts/infra_compendium/load.py`. Rebuild: `scripts/infra_compendium/run_all.sh`. Env: `/home/plm_env/bin/python`.
- **Inputs**:
  - `/workspace/data/01–03` (compendium: 1,937 conditions, 178 PMIDs, 189,962 sites, 17.0 M log2FC values)
  - PSP 2025-05 files (Kinase_Substrate_Dataset, Phosphorylation_site_dataset, Regulatory_sites)
  - Klaeger 2017 Table S2 (high-confidence kinobeads targets, Kd_app)
  - `kinhib_inhib_comb__1.csv` (Hijazi 2020)
- **Read-only resources reused**:
  - kinaseActivityMain `mvp/` (core-v3 code and model artefacts, copied to `scripts/infra_compendium/kinact_mvp_pkg/`)
  - `cache/libraries/curated_union.parquet` (the evaluated curated library)
  - `results/ground_truth/{drug_targets_manual,drug_aliases,condition_doses,stimulus_kinases}.tsv` (curated drug→target map, aliases, doses, stimulus kinases)
  - `external/uniprot` annotations (gene symbols)
  - kinaseSubMain `results/labels/kinase_list.tsv` (Manning/KinHub families) and `results/release/{site_top5,kinase_top200}.tsv.gz` (predicted substrates; CC BY-NC-SA)
- **Downloaded by me**: IDG understudied-protein list, `github.com/druggablegenome/IDGTargets` (`IDG_TargetList_CurrentVersion.json` sha256 `0bb0c011…`, `IDG_List_20210120.xlsx` `cfcb69a6…`), saved to `/workspace/cache/compendium/external/idg/` on 2026-10-29. It has 134 kinases; 126 map to the 503-kinase universe.
- **Kinase activity** (kinase × condition):
  - `corev3`: z, z_core, z_sig, n_substrates. Edges `kinact_mvp_pkg/model/edges.tsv.gz` (sha256 `d079a72b…`) = curated_union + r2np top-100. This is the library of the evaluated core-v3. The MVP `build_edges.py` rebuild (PSP + OmniPath) is kept at `cache/compendium/edges_mvp_rebuild.tsv.gz` and shares 91% of the curated edges.
  - `psp_z`: PSP-only Stouffer z of robust-standardised substrates, at least 5 substrates. This equals CPTAC's "meanz" × √n, so the two rank kinases differently when their substrate counts differ.
- **Licences**: PSP-derived tables (edges, `psp_*` columns, predictions) are CC BY-NC-SA 3.0, non-commercial. Klaeger and Hijazi supplements are for citation only; do not redistribute.
