# val_E: Finding E and I24 offsets re-tested with a cold-ischaemia covariate

Scripts: `/workspace/scripts/val_E/v01_isch.py`, `v02_main.py`, `v03_offsets.py`. Tables are in this directory. Per-tumour score: `ischaemia_scores.tsv`.

## Bottom line
Finding E is robust to the ischaemia covariate. So are the KRAS offset and the LUAD EGFR-mutant result. The ischaemia score is essentially uncorrelated with receptor protein and with KRAS status in the cohorts that carry these analyses. It therefore cannot confound them, and adding it moves the estimates by less than their CIs.

## Method
- **Score.** The audit_I10 per-tumour score was not saved, so I recomputed it with that audit's own code (`audit_I10/common.py`: `ischaemia_set`, `cohort`, `standardise`, `score`). It is the signed mean z of the top-200 Mertins 2014 cold-ischaemia sites, on `input_tc_padj`, tumours only.
- **Circularity fix.** The 200-site set contains ERK1 pT202 (P27361_T202), which is part of the ERK outcome. The primary score drops all sites on EGFR, ERBB2, MET, SHC1, GAB1, PTPN11 and ERK1/2. It correlates r ≥ 0.997 with the raw score in every cohort.
- **E recompute.** I followed the audit_I24 code (`build`, `a03` joint model, `decomp`). The `tc` variant is the one carrying the −0.15. The model is ERK pT/pY ~ RTK pY + RTK protein + purity + E2F + isch, as OLS on z-scores with a DerSimonian–Laird random-effects meta, 7 cohorts and n = 613. The baseline reproduces the audit exactly (−0.151).
- **Offsets.** These use the b_I24 T4 'main' model (Huber, RTK>ERK, `tc_padj`, mutation lesion) with isch added as an additive covariate. The baseline reproduces +0.674.
- **Tissue source site is not available.** The CPTAC clinical and biospecimen tables have no TSS or ischaemia-time column. The only proxy is the numeric case-ID prefix, which exists for BRCA, COAD and OV only. Only BRCA and COAD enter E. As fixed effects they are uninformative (2 cohorts, CI −0.55 to 0.27). I report them for completeness only.

## Task 1: before and after

### Claim E, abundance coefficient (ERK ~ RTK protein at fixed RTK pY)
| specification | meta coef [95% CI] | p | cohorts negative | I² |
|---|---|---|---|---|
| Baseline (I24 audit) | −0.151 [−0.248, −0.053] | 0.002 | 6/7 | 0.05 |
| **+ isch (primary)** | **−0.158 [−0.254, −0.062]** | **0.001** | 6/7 | 0.07 |
| + raw isch (includes ERK1 T202) | −0.158 [−0.252, −0.063] | 0.001 | 6/7 | 0.05 |
| + isch computed from `tc` input | −0.154 [−0.247, −0.062] | 0.001 | 6/7 | 0.00 |
| + isch + isch² | −0.152 [−0.246, −0.058] | 0.002 | 6/7 | 0.00 |
| + isch + site-ID proxy (BRCA, COAD only) | −0.141 [−0.553, 0.271] | 0.50 | 1/2 | – |

Per cohort (base → +isch): BRCA 0.035 → 0.026, CCRCC −0.116 → −0.134, COAD −0.160 → −0.130, GBM −0.089 → −0.112, HNSCC −0.178 → −0.148, LSCC −0.125 → −0.149, LUAD −0.377 → −0.381. The LUAD estimate is the largest and was untouched. BRCA is the one positive cohort and remains so.

The coefficient on RTK pY goes from 0.116 [−0.002, 0.234] to 0.133 [0.017, 0.248].

**Negative control.** To check for over-adjustment or artefact, I adjusted for 30 random 200-site scores (random signs, pathway proteins excluded) in place of the ischaemia score. The coefficient is −0.146 ± 0.018 (range −0.172 to −0.103). The ischaemia score changes it by no more than a random score does, in the direction of strengthening it.

### Second half of E: "once abundance is removed, RTK-driven adaptor pY transmits about as well as RTK-independent pY"
IV decomposition (500-bootstrap robust SEs, 7 cohorts):

| variant | spec | RTK-driven (IV) | RTK-independent | diff | diff p |
|---|---|---|---|---|---|
| `tc` | base | 0.14 [−0.07, 0.36] | 0.45 [0.32, 0.57] | −0.32 | 0.006 |
| `tc` | +isch | 0.17 [−0.02, 0.36] | 0.44 [0.30, 0.57] | −0.29 | 0.006 |
| `tc_padj` | base | 0.31 [0.13, 0.49] | 0.44 [0.32, 0.56] | −0.13 | 0.18 |
| `tc_padj` | +isch | 0.35 [0.18, 0.51] | 0.42 [0.28, 0.55] | −0.08 | 0.35 |

The abundance-adjusted (`tc_padj`) transmission gap shrinks slightly with isch (0.13 to 0.08). The gap is still not significant, as in the audit. The `tc` gap, which is the abundance part, is unchanged.

### KRAS offset (RTK→ERK, `tc_padj`, 3 cohorts with ≥5 mutants, 51 mutants)
| spec | meta offset [95% CI] | p |
|---|---|---|
| Base (I24) | +0.674 [0.337, 1.012] | 9e-5 |
| + isch | +0.598 [0.239, 0.958] | 0.001 |

By cohort (base → +isch): LUAD 0.78 → 0.73, COAD 0.64 → 0.55, LSCC 0.04 → −0.07.

### LUAD EGFR-mutant (quick check)
- **Offset at fixed RTK input** (b_I24 model; EGFR-mut n = 39 vs all non-mut): −0.95 ± 0.18 → −0.99 ± 0.16.
- **ERK level, EGFR-mut (n = 36) vs driver-negative (n = 27)**, Huber, adjusted for purity and E2F. "Driver-negative" means no mutation in KRAS, EGFR, BRAF, NRAS, HRAS, ERBB2, MET, RIT1 or MAP2K1. EGFR-mut with a co-mutation was dropped. This definition is mine; the I24 tables contain no pre-existing driver-negative contrast. Result: −0.76 SD ± 0.25 (z ≈ −3.1) → −0.77 ± 0.24 with isch.

## Task 2: does the ischaemia score track receptor protein or KRAS?
Spearman ρ with the score (`isch_correlations.tsv`):

| variable | median \|ρ\| over cohorts | range |
|---|---|---|
| EGFR protein | 0.14 | −0.26 to 0.23 |
| ERBB2 protein | 0.06 | −0.16 to 0.17 |
| RTK protein composite (the covariate in E) | 0.08 | −0.19 to 0.17 |
| RTK pY | 0.06 | −0.15 to 0.06 |
| Adaptor pY | 0.12 | −0.12 to 0.24 |
| ERK pT/pY (outcome) | 0.33 | 0.05 to 0.43 |

- **KRAS.** In the three cohorts used for the offset the correlation is ρ = −0.13 (COAD), 0.01 (LSCC) and 0.08 (LUAD). The KRAS-mutant minus wild-type difference in score is −0.07, 0.14 and 0.18 SD (Mann–Whitney p = 0.22, 0.92 and 0.39).
- **Outside the analyses, UCEC** shows KRAS-mutant tumours with +0.73 SD higher score (p = 9e-4). UCEC has no RTK stage, so it does not enter E or the offset.
- **Why it does not confound.** Only the outcome (ERK pT/pY) is substantially tied to the score. The score has no consistent link to the receptor protein or to KRAS in the contributing cohorts. A confounder needs to relate to both exposure and outcome. The tie to the outcome is the audit_I10 finding, and it matters for levels of MAPK output. It does not transmit into the adjusted contrasts here.

## Task 3: verdicts
| claim | before | after | verdict |
|---|---|---|---|
| E: at fixed RTK pY, more receptor protein goes with lower ERK pT/pY (`tc`) | −0.15 [−0.25, −0.05], 6/7 | −0.16 [−0.25, −0.06], 6/7 | **Robust** |
| E: once abundance is removed, RTK-driven adaptor pY transmits about as well as RTK-independent pY (`tc_padj`) | 0.31 vs 0.44, p = 0.18 | 0.35 vs 0.42, p = 0.35 | **Robust** (still not significantly different; if anything closer) |
| KRAS mutation raises ERK at a given RTK input (+0.67 SD) | +0.67 [0.34, 1.01] | +0.60 [0.24, 0.96] | **Robust** (about 11% smaller, CIs overlap, 3 cohorts, LSCC has only 5 mutants) |
| LUAD EGFR-mut ERK below driver-negative | −0.76 ± 0.25 | −0.77 ± 0.24 | **Robust** (small n: 36 vs 27) |

## Caveats
- **Proxy score.** It is an expression-style footprint from three tissue experiments, not a measured ischaemia time. Measurement error in it would leave residual confounding. The negative result is more convincing for the lack of any correlation with exposures (Task 2) than for the unchanged coefficients alone.
- **Over-adjustment.** The score also tracks real MAPK output, since it correlates 0.33 with ERK. Adjusting for it can therefore over-correct, and that would push the coefficients toward 0. They did not move, so this does not weaken the conclusion.
- **Real TSS or cold-ischaemia time** (from the CPTAC biospecimen metadata) is not in the assets. If it can be obtained, this check could be repeated.
- **KRAS offset.** It rests on 3 cohorts.
