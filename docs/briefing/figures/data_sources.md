# Data sources for data figures (fig3–fig7)

All panels re-plot saved result tables. No analysis was re-run. Scripts: `docs/briefing/figures/src_data/` (`style.py` holds the shared house style; one script per figure). Run with `/home/plm_env/bin/python <script>` from that directory, `OMP_NUM_THREADS=4`.

The only derived file is `src_data/derived_menyhart2023_psp_ck2_sites.tsv`. It is written by `s00_menyhart_psp_ck2.py`, which reads the downloaded Menyhart 2023 supplementary tables (`cache/audit_F/dl/menyhart2023/mmc4.xlsx`, `mmc5.xlsx`; run with `python -I`) and keeps the PhosphoSitePlus-annotated CK2 substrates. It repeats audit_F's own sanity check (`scripts/audit_F/t2_query_menyhart.py`) and saves it per site. Its medians match `results/audit_F/AUDIT.md`: M/L −0.57 to −1.33; SGC-CK2-1 H/M +0.78 / +1.16.

## fig3_kpna4_ck2 — K1: KPNA4 S60 is a CK2 site
| Panel | What is shown | Source file(s) | Columns / filter |
|---|---|---|---|
| a, CX-4945 points | KPNA4 S60 robust z, rank and percentile in each CX-4945 experiment (HeLa, NTERA2, HL60; MCF7 not measured) | `results/audit_D/tables/kpna4_ck2i_rank.tsv` | `agent == CX4945`, `z`, `rank`, `n_sites`, `pct` |
| a, ERK points | KPNA4 S60 z in the 25 units with verified ERK-pathway inhibition; median and count below −2 | `results/val_I03/tables/t3_mapk_concordance.tsv` | `pathway_down == True`, `KPNA4_S60` |
| b, KPNA4 S60 and fast-CK2-set sites | log2 M/L (WT CK2α + inhibitor vs DMSO), H/L (resistant CK2α + inhibitor vs DMSO), H/M rescue, FDR | `results/audit_F/task2_menyhart2023_module_sites.tsv` | quantified rows, `multiplicity == 1` |
| b, known CK2 substrates | median M/L and H/L of PSP-annotated CK2 substrates in each arm | `docs/briefing/figures/src_data/derived_menyhart2023_psp_ck2_sites.tsv` (from `cache/audit_F/dl/menyhart2023/mmc4.xlsx`, `mmc5.xlsx`) | median per inhibitor × time |
| c | Sequence window HEDICED**S**DIDGDYR | `results/audit_D/tables/cx4945_specific_module.tsv` (`window` column for O00629_S60); motif text from `results/audit_D/AUDIT.md` | — |

## fig4_tissue_wiring — K2: perturbation programmes in human tissue
| Panel | What is shown | Source file(s) | Columns / filter |
|---|---|---|---|
| a, lines | Coherence (mean signed site-to-site r) of MEKi, EGF, TNF and DNA-damage footprints in normal adjacent tissue (NAT) and tumour, 6 paired cohorts | `results/b_Q8/tables/T1_coherence.tsv` | `cond == padj_pc3`, `design == paired`; `coh_N`, `coh_T` |
| a, grey diamonds | Matched random-site null coherence (median over cohorts) | same | `null_N`, `null_T` |
| a, z labels, meta z | Median z vs matched random sites, NAT → tumour; DNA-damage meta z | `results/b_Q8/tables/T2_coherence_summary.tsv` | `cond == padj_pc3`; `z_N_med`, `z_T_med`, `meta_z_perm` |
| b | Sign agreement with the compendium for the 529 compendium relationships, NAT vs tumour, per cohort (protein-adjusted) | `results/b_Q8/tables/T7_edge_global_preservation.tsv`; 529 = `label == edge` count in `T3_compendium_edges.tsv` | `cond == padj`, `label == edge`; `signconc_N`, `signconc_T`; median over the 6 paired cohorts |
| c | Tumour coherence with no adjustment, after adjusting for measured ischaemia time (t + log t), and after adjusting for the proxy ischaemia score; 7 cohorts with GDC times | `results/val_B/T3_coherence.tsv` | `adj` in {`none`, `measured t + log t`, `footprint ISCH score`} |

## fig5_portability — K3: primary vs reverse responses
| Panel | What is shown | Source file(s) | Columns / filter |
|---|---|---|---|
| a | Cross-study replication log odds ratio (95% CI) by effect-size bin, decreases vs increases, inhibitors vs stimulations | `results/audit_I19/T_fig_mirror.tsv` | all rows; study counts (41, 52) from `results/audit_I19/AUDIT.md` |
| b | Leave-one-study-out footprint AUROC using decrease sites vs increase sites (auditor's own code), per footprint | `results/audit_I10/A3_loso_own.tsv` | rows with both `down` and `up`; `kind` splits inhibitor vs stimulus/stress |

Note on (b): `results/b_I22/tables/` has no per-site table for the MCF7 Hijazi vs Ochoa contrast; the numbers exist only as text in `b_I22/NOTES.md` and `SUMMARY.md`. As the brief allowed, panel b therefore uses the footprint-transfer fallback from audit_I10.

## fig6_cascade — K4: offsets not gains; buffered receptor abundance
| Panel | What is shown | Source file(s) | Columns / filter |
|---|---|---|---|
| a | Schematic (no data) | — | — |
| b | Random-effects meta gain per cascade link (standardised slope, 95% CI, number of cohorts). The "expected" tick is the product of the RTK→adaptor and adaptor→ERK gains | `results/b_I24/tables/T2_gain_meta.tsv` | `variant == tc_padj`; `std_slope_meta`, `std_lo`, `std_hi`, `k` |
| c | Gain change (`int_meta`) and offset change (`lev_meta`) with 95% CI and BH q for selected pre-declared lesion tests; count of significant gain tests | `results/b_I24/tables/T4_modifiers_meta.tsv` | `model == main`, `primary == True` (162 tests); `int_*`, `lev_*`, `*_q_primary` |
| c, shaded band | Smallest gain change detected with ≥80% power in spike-in (0.15: DUSP6 94%, PTEN 84%) | `results/b_I24/tables/T6_spikein_power.tsv` | `power >= 0.8` |
| d | ERK pT/pY coefficient on receptor protein at fixed receptor pY, per cohort (±1.96 SE) and DerSimonian–Laird meta | `results/val_E/E_joint_by_cohort.tsv`, `results/val_E/E_joint_meta.tsv` | `spec == "base (purity, prolif)"`; `b_prot`, `se_prot`; meta `prot`, `prot_lo`, `prot_hi`, `prot_p`. Matches audit_I24 (−0.15 [−0.25, −0.05], p = 0.002, 6/7 cohorts negative; `results/audit_I24/AUDIT.md` §2) |

## fig7_caveated — three caveated observations
| Panel | What is shown | Source file(s) | Columns / filter |
|---|---|---|---|
| a | Class-median z under catalytic mTOR, AKT and S6K inhibitors for candidates (LAMTOR1 S98, S97; RRAGC S15), direct mTORC1 controls and S6K controls; number of mTOR-inhibitor units; LAMTOR1 S98 percentile | `results/audit_D/tables/class_effects_loso.tsv` | `cls` in {`cat_mTORi`, `AKTi`, `S6Ki`}; `med_z`, `n_units`, `pct_vs_sites` |
| b | Split-by-lab Spearman ρ of response speed for substrates vs random sites vs noise-matched random sites, per kinase family (95% CI) | `results/audit_H/h1_lab_splits.tsv` (substrate, random), `results/audit_H/h1_noise_matched_random.tsv` (noise-matched random) | `analysis == lab_split`; `set == random_noise_matched`, `level == lab` |
| c | Target-substrate potency shift (pEC50 − Davis 2011 ATP-free pKd) vs the Cheng–Prusoff band for 1–5 mM ATP; implied [ATP] per drug | `results/b_I17/atp_drug_table.tsv`, `results/b_I17/atp_implied.tsv` | rows with `shift_tsub_davis`; band `CP_pred_5mM`..`CP_pred_1mM`; `implied_ATP_mM` |
