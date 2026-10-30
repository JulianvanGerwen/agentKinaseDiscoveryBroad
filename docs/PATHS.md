# Paths and placeholders

The repository was cleaned of machine-specific paths before release. Some scripts still reference the original layout, so a reader must set a few paths before re-running them.

## Placeholders

- `<HACKATHON_ROOT>/` stands for the directory that holds the two sibling tool repositories `kinaseSubMain/` and `kinaseActivityMain/`. Scripts read inputs from `<HACKATHON_ROOT>/kinaseSubMain/external/...` (CPTAC PDC files, DepMap, GTEx, HPA, HGNC, GO, UniProt, kinase library, interactomes, Ochoa 2020, KinCore, Davis 2011), from `kinaseSubMain/results/{release,labels}/`, and import `kinaseActivityMain/mvp`. Replace the placeholder with your own checkout location (for example by a global find and replace, or a symlink).
- `<SCRATCH>` stands for a temporary working directory (used only by `scripts/val_1/v2.py` and `v2b.py` for a pickle hand-off). Use any writable directory.

## Paths a reader must provide

- `/workspace/` is the project root assumed by the scripts. Clone to `/workspace`, or symlink it, or find and replace.
- `/workspace/data/` holds the provided phosphoproteomics compendium (`01` to `03` files) and PhosphoSitePlus files (`Kinase_Substrate_Dataset`, `Phosphorylation_site_dataset`, `Regulatory_sites`) and `Klaeger/` tables. Not included (licence-restricted); obtain them yourself.
- `/workspace/cache/` and `results/*/cache/` hold downloads and intermediates (parquet, npy and so on). Not included; scripts recreate them. `cache/**/dl/` folders hold downloads from GDC, GDSC, DepMap and CPTAC/PDC (see `NOTICE.md`).
- `/home/plm_env/bin/python` is the Python environment used in `run_all.sh` scripts (Python 3.12, pandas, pyarrow, torch). It is a generic environment path; point it at your own interpreter. It occurs in: `DATA.md`, `RESEARCH_TREE.md`, `results/infra_compendium/SUMMARY.md`, `results/infra_cptac/SUMMARY.md`, `scripts/b_I03/run_all.sh`, `scripts/b_I26/run_all.sh`, `scripts/b_I28/run_all.sh`, `scripts/infra_compendium/run_all.sh` and `scripts/infra_compendium/logs/01_build_matrix.log`.
- `/home/.cache/gseapy/Enrichr.GO_Biological_Process_2023.gmt` is read by `scripts/b_I10/06_module_annot.py` (gseapy cache; download the Enrichr GO Biological Process 2023 library).

## Environment variables

Several scripts take optional settings via environment variables, all with defaults: `TAG`, `SET_TAG`, `MODE` (`lolo` default), `LABELS`, `PCS`, `NPERM`, `MARG`, `TOP`, `ONLY`, `ICA_CENTER`. See the head of each script.

## Excluded from the repository

Patient-level tables (CPTAC/GDC case and sample identifiers) were removed: `results/b_I29/tables/T7_cptac_conf_tumour_table.tsv`, `results/val_B/gdc_ischaemia_all.tsv`, `results/val_E/ischaemia_scores.tsv`. Scripts that write them will regenerate them from the source data. Large tables and all binary matrices are excluded as listed in `.gitignore`.
