# Phosphoproteomics compendium: analysis-ready assets (infra_compendium, 2026-10-29)

**Rebuild:** `scripts/infra_compendium/run_all.sh` (about 2 min, under 8 GB RAM, 8 threads). **Env:** `/home/plm_env/bin/python` (py3.12, pandas 3.0, pyarrow 25, torch 2.6). **Loader:** `scripts/infra_compendium/load.py`.

## Files (`/workspace/results/infra_compendium/`)
| file | content / schema |
|---|---|
| `site_matrix.parquet` | 189,962 sites × 1,937 conditions. Column `site` (`{UniProt}_{Res}{Pos}`) + one float32 log2FC column per `conditionID`. NaN = not measured. Long cache: `cache/compendium/quant_long.parquet` |
| `sites.tsv` | site, uniprot, residue, position, gene (99.4% mapped, UniProt 2026_03), window, n_conditions, n_studies, `psp_kinases_human`, `psp_regulatory`, `pred_top1_kinase`/`_assignment_prob`/`pred_top3_kinases` (kinaseSubMain release) |
| `condition_coverage.tsv` | per condition: n_sites, n_psp_ks_sites, n_proteins, median, mean, sd, mad, max_abs, frac_abs_gt1, duplicate_of |
| `conditions_annotated.tsv` | conditionID, study_pmid, perturbation_class/subclass, agents, kinase_inhibitors, other_drugs, control_agents, background_agents, contrast_reversed, dose_nM, inhibitor_primary_targets(_in_universe), inhibitor_profiled_targets, target_confidence, genetic_target_kinases (`GENE(+/-)`), stimulus, stimulus_expected_up_kinases, cell_type, cell_line_norm, sample_type, time_min, time_source, time_series_id, time_rank, dose_series_id, qc_flag |
| `drug_targets.tsv` | 209 agents (154 kinase inhibitors): curated primary targets, Klaeger Kd_app ≤100 nM / ≤1 µM (high-confidence targets, Table S2), Hijazi targets, target_confidence |
| `activity_long.parquet` | condition, kinase, method (`corev3` / `psp_z`), z, n_substrates, mean_log2fc (psp_z), z_core, z_sig, mean, sd (core-v3) |
| `activity_{corev3,pspz}_{z,nsub}_wide.parquet` | kinase × condition (float32) |
| `kinase_annotation.tsv` | 503 kinases: Manning group/family/subfamily (KinHub), pseudokinase, `idg_dark_kinase`, PSP substrate counts (total / in compendium), core-v3 edge counts, per-method coverage |
| `sanity_inhibitor_targets{,_summary}.tsv` | inhibitor-target sanity check |

## Key counts
- 178 studies; 41% of conditions come from decryptM (PMID 36926954, 790) and 9% from Hijazi (31959955, 183). The matrix is 4.6% filled. A median site is measured in 11 conditions; 32,758 sites in ≥100 conditions; 3,724 in ≥50%. Sites per condition: median 8,560 (1,981–29,577).
- **Conditions by class (n conditions / n studies):**
  - kinase_inhibitor: 850 / 55 (154 agents; target confidence high 522, medium 167, low 161)
  - other: 564 / 63 (306 are rituximab dose × time in decryptM, the rest drugs, infection, cell cycle and GPCR ligands)
  - stress: 142 / 27
  - cytokine: 97 / 11
  - time_course (later vs earlier time of the same treatment): 92 / 3
  - disease_vs_control: 75 / 21
  - genetic: 69 / 31
  - growth_factor: 48 / 13
- Time is parsed for 1,664 conditions. 714 conditions fall into 172 time series and 835 into 86 dose series. qc_flag: 120 low_coverage, 12 scale_suspect, 7 duplicate_vector.
- **Activity, core-v3:** the evaluated library (curated_union + r2np top-100; 27,598 curated + 46,866 predicted edges). All 503 kinases get a score in every condition, through family borrowing. 352 kinases have ≥5 measured substrates in ≥50% of conditions.
- **Activity, psp_z:** PSP-only, robust per-condition z, Stouffer score, ≥5 substrates. 189 kinases are scored at least once, 66 in ≥50% of conditions. Correlation with core-v3: r = 0.58.
- **Sanity check:** 488 inhibitor conditions from 39 studies (high/medium confidence, qc ok, dose ≥100 nM). Median AUROC of targets being decreased:
  - core-v3: 0.95 (targets z<0: 86%)
  - core-v3 z_core only: 0.85 (79%)
  - psp_z: 0.86 (78%; 354 conditions)

## Caveats
- **Batch.** Study equals batch. Use `study_pmid` for grouping and CV, and never pool z across studies without it. decryptM dominates by condition count.
- **Contrast is condition vs its `Control` column.** It may be reversed (`contrast_reversed`; e.g. "DMSO vs ceritinib"), a comparison of two inhibitors (`inhibitor_vs_other_inhibitor`), an effect on top of a background drug (`background_agents`, e.g. TNF under an inhibitor), or a comparison with an earlier time point (`time_course`).
- **Unscored kinases.** For core-v3, z with `n_substrates`=0 comes only from family borrowing and the prior. Filter on n_substrates ≥1 (or ≥5).
- **In-sample signatures.** The core-v3 signature term (`z_sig`, 55 kinases) and its calibration were fitted on development ground-truth conditions *of this compendium*. Its sanity numbers are therefore optimistic. Use `z_core`, or psp_z, for independent tests.
- **Leakage filter.** The library drops edges that cite the benchmark studies, which makes it conservative.
- **Agent parsing.** Agents are parsed with rules and aliases. Unknown compounds are kept as text (AZ123, TRC-794, "SM", EnkaninE). Times in qPTM rows are often missing. Doses are mostly defaults (10 µM) outside decryptM and Hijazi.

## Loading
```python
import sys; sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
C = load.conditions(); X = load.sites_matrix(C.conditionID[C.study_pmid == "31959955"])
Z = load.activity("corev3"); N = load.activity("corev3", "nsub"); Z = Z.where(N >= 1)
P = load.activity("psp_z"); K = load.kinases(); D = load.drug_targets()
```
