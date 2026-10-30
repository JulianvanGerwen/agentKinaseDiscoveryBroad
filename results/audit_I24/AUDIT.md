# Audit of I24 candidate E: "non-transitive RTK → adaptor → ERK transmission"

Auditor scripts: `/workspace/scripts/audit_I24/` (common_a.py, a01–a05). Tables: `/workspace/results/audit_I24/tables/` (A0–A4b). Logs: `a01–a04.log`. Figure: `figures/AF1_decomposition_specifications.png`. I did not modify any b_I24 file.

## Verdict: QUALIFIED (fails as worded)

The 2SLS gap reproduces in unadjusted pY data (`tc`). It is not a general non-transitivity of the cascade, for three reasons:
- It is not significant on protein-adjusted (stoichiometry) data.
- It disappears in unadjusted data once RTK protein abundance is adjusted for.
- In stoichiometry units, RTK→ERK matches the product of the two edges.

What survives is narrower and largely known. Tumour-to-tumour differences in RTK pY that come from **receptor abundance** (amplification or overexpression) raise adaptor pY but do not reach ERK. Regression dilution is not the explanation, and the direction of that bias was checked rigorously (§3).

**Defensible wording:** "In CPTAC tumours, EGFR/ERBB2 pY that tracks receptor protein abundance does not propagate to ERK1/2 activation-loop phosphorylation. At fixed receptor pY, higher receptor protein goes with lower ERK pT/pY (−0.15 SD [−0.25, −0.05], p = 0.002, 6/7 cohorts). Once receptor abundance is removed, RTK-instrumented adaptor pY transmits to ERK about as well as RTK-independent adaptor pY: 0.31 [0.13, 0.49] vs 0.44 [0.32, 0.56], p = 0.16. Acute EGF/EGFR-inhibitor perturbations in cell lines are transitive. This fits the known decoupling of chronic receptor overexpression from steady-state MAPK output (robustness/feedback). It is not evidence of a general break in the cascade."

**Most convincing figure:** `figures/AF1_decomposition_specifications.png`. It is a forest plot of the RTK-driven part vs the RTK-independent part across six specifications plus the cell-line contrast. The gap is present only in the two rows that keep RTK abundance in the instrument.

## 0. Problems with how the claim is built (from I24's own tables)
- **Mixed variants.** The headline gains (RTK→adaptor 0.41, adaptor→ERK 0.41, RTK→ERK 0.10) come from `tc_padj`. The decomposition with p = 4e-4 comes from `tc`. On `tc_padj`, I24's own decomposition gives 0.28 vs 0.44, p = 0.07.
- **The transitivity test with the quoted gains does not reject transitivity.** Expected RTK→ERK = 0.424 × 0.393 = 0.17, and the observed 0.10 [0.02, 0.18] contains it. Non-transitivity holds only in `tc`: 0.357 × 0.386 = 0.14 expected vs 0.03 [−0.05, 0.10] observed.
- **"RTK-independent" means independent of 1–6 measured EGFR/ERBB2 sites.** Other RTKs, SFKs and global pY are all inside the residual (§4).
- **"10 cohorts" is really 7.** The RTK stage is missing in OV, PDAC and UCEC. In BRCA it is a single site (ERBB2 Y1248, measured in 60% of tumours); in COAD it is a single site (EGFR Y1092, 44%).
- **The "2SLS" is a Hausman-type contrast.** The coefficient on `a_rtk` is reduced form / first stage. The coefficient on `a_oth` is the ordinary partial regression slope of ERK on adaptor pY controlling for RTK.

## 1. Independent recompute (own code, same site definitions; OLS standardized, 1000-bootstrap robust SE, DerSimonian–Laird meta)

| | k | first-stage F (median, min) | RTK-driven part (IV) | RTK-independent part | diff p |
|---|---|---|---|---|---|
| I24, `tc` | 7 | – | 0.12 [−0.05, 0.30] | 0.46 [0.34, 0.57] | 4e-4 |
| audit, `tc` | 7 | 13.3, 5.9 | 0.13 [−0.08, 0.34] | 0.45 [0.32, 0.57] | 0.004 |
| I24, `tc_padj` | 7 | – | 0.28 [0.12, 0.45] | 0.44 [0.32, 0.55] | 0.07 |
| audit, `tc_padj` | 7 | 18.5, 5.2 | 0.31 [0.13, 0.49] | 0.44 [0.32, 0.56] | 0.16 |

Per-cohort IV and RTK-independent estimates, `tc`:

| cohort | IV | RTK-independent |
|---|---|---|
| BRCA | −0.18 | 0.33 |
| CCRCC | −0.11 | 0.63 |
| COAD | 0.48 | 0.54 |
| GBM | −0.15 | 0.27 |
| HNSCC | 0.00 | 0.29 |
| LSCC | 0.11 | 0.62 |
| LUAD | 0.38 | 0.46 |

Only LSCC is individually significant (p = 0.004). COAD and LUAD show no gap.

The audit gains are close to I24's:

| edge | `tc_padj` | `tc` |
|---|---|---|
| RTK→adaptor | 0.42 | 0.36 |
| adaptor→ERK | 0.39 | 0.40 |
| RTK→ERK | 0.10 | 0.03 |

The point estimates reproduce. The p-value is about 10× weaker with bootstrap SEs: the RLM model SE understates the heavy-tailed sampling error of the ratio.

The result depends on the outcome isoform:
- ERK2-only outcome (`tc`): 0.25 vs 0.40, p = 0.10.
- ERK1-only outcome (`tc`): 0.04 vs 0.41, p = 0.001.
- Within the adaptor stage, the gap on `tc_padj` is null for SHC1 alone (0.33 vs 0.38, p = 0.61) and for PTPN11/GAB1 alone (0.36 vs 0.35, p = 0.83).

AKT substrates (my recompute): `tc` 0.24 vs 0.40 (p = 0.053); `tc_padj` 0.27 vs 0.41 (p = 0.35).

## 2. Instrument strength and composition; amplified receptors
- **Composition of the "RTK pY" instrument:**
  - EGFR pY only: CCRCC, COAD, GBM, LSCC, LUAD.
  - ERBB2 Y1248 only: BRCA.
  - EGFR + ERBB2: HNSCC.
  - MET is never measured.
  (`A0_stage_sites.tsv`)
- **First-stage F (`tc`):**

  | BRCA | CCRCC | COAD | GBM | HNSCC | LSCC | LUAD |
  |---|---|---|---|---|---|---|
  | 10.2 | 5.9 | 11.8 | 13.3 | 13.4 | 22.7 | 18.0 |

  CCRCC is weak. On `tc_padj`, BRCA (5.2), COAD (9.4) and CCRCC (9.7) are below 10. Weak-instrument bias in 2SLS points toward OLS, so it would shrink the gap, not create it. Restricting to cohorts with F ≥ 10: `tc` 0.15 vs 0.41, p = 0.009; `tc_padj` 0.31 vs 0.40, p = 0.32.
- **Amplification drives the instrument where it exists:**
  - BRCA: excluding the 16 ERBB2-amplified tumours collapses the first stage (F = 0.16 on `tc`, 0.001 on `tc_padj`). The BRCA instrument is entirely ERBB2 amplification.
  - GBM: excluding the 47 EGFR-amplified tumours drops F to 2.2.
  - After excluding amplified tumours, only LSCC, LUAD (and HNSCC on `tc_padj`) have F ≥ 10: `tc` 0.17 vs 0.60 (p = 0.075); `tc_padj` 0.36 vs 0.44 (p = 0.46).
- **Per receptor, EGFR pY only:** `tc` 0.16 vs 0.46 (p = 0.008); `tc_padj` 0.32 vs 0.46 (p = 0.21).
- **Phosphosite level tracks receptor protein:** in unadjusted data EGFR pY vs EGFR protein has Spearman ρ = 0.80 (GBM), 0.69 (HNSCC), 0.64 (LSCC), 0.65 (LUAD).
- **Decisive dissection (`A2_protein_dissection_meta.tsv`, `tc`):**
  - Adding adaptor + ERK protein as covariates leaves the gap intact: 0.12 vs 0.44, p = 0.002.
  - Adding RTK protein alone removes it: 0.29 [0.05, 0.52] vs 0.44, p = 0.20.
  - Splitting RTK pY into an abundance-predicted component and a stoichiometry component (`A3`): the stoichiometry instrument gives 0.28 [0.05, 0.52] vs 0.41 (p = 0.19). The abundance component is too weak to estimate (F 0.3–6).
  - Joint reduced form: ERK ~ RTK pY + RTK protein gives RTK protein −0.15 [−0.25, −0.05] (p = 0.002; negative in 6/7 cohorts) and RTK pY +0.12. The analogous adaptor model gives RTK protein −0.10 (p = 0.06).

So the non-transitivity is the abundance part of RTK pY.

## 3. Measurement / regression dilution (main alternative) — it does not explain the gap; its bias runs the other way
- **Analytic argument.** Classical noise in the instrument (RTK pY) shrinks the first stage and the reduced form by the same factor, so the IV ratio stays unbiased; noise only widens the CI. Noise in the mediator (adaptor pY) attenuates the RTK-independent coefficient, the OLS partial slope. Classical noise therefore biases the contrast *in favour of* the RTK-driven part.
- **Split-half reliabilities (Spearman–Brown, adjusted for covariates; `A2_split_half_reliability.tsv`, `tc`):**

  | cohort | RTK | adaptor | ERK |
  |---|---|---|---|
  | GBM | 0.90 | 0.67 | 0.71 |
  | HNSCC | 0.92 | 0.32 | 0.54 |
  | LSCC | 0.96 | 0.76 | 0.79 |
  | LUAD | 0.94 | 0.71 | 0.90 |
  | CCRCC | 0.74 | 0.85 | 0.86 |

  For BRCA and COAD the RTK stage is a single site, so its reliability cannot be split. On these numbers RTK pY is the *most* reliable stage, not the noisiest.
- **Reliability-matched simulation null (`A2_sim_summary.tsv`).** I simulated full transitivity using each cohort's n, first stage, total gain and the observed reliabilities, 300 replicates per cohort. The expected IV − residual difference is **+0.52** (`tc`; 95% range 0.24–1.52) and +0.46 (`tc_padj`). The observed values are −0.37 and −0.21, with P(null ≤ observed) = 0.010 and 0.007. Classical noise predicts the opposite sign.
- **Reliability-matched comparator: an IV within the adaptor stage.** Instrumenting PTPN11/GAB1 pY with SHC1 pY, or the reverse, gives adaptor→ERK = 0.69 [0.52, 0.87] and 0.69 [0.54, 0.84] in the 7 cohorts. So adaptor variation identified by an instrument with similar reliability transmits strongly; only the RTK-pY-identified variation (`tc`) does not. Caveat: a within-adaptor IV is inflated by any factor the halves share (global pY, SFK), so treat it as an upper bound.
- **Correlated (non-classical) error is the remaining measurement route.** I proxied it with a global pY factor (`bgPY`, the mean of about 60–210 non-pathway pY sites). The residual adaptor pY correlates with bgPY (ρ 0.1–0.6). Adjusting for bgPY lowers the first stage (0.35 → 0.25) and sets RTK→ERK to about 0, but the gap persists in `tc` (−0.08 vs 0.35, p = 0.009) and is not significant on `tc_padj` (p = 0.10). Part of the small `tc_padj` RTK→ERK gain (0.10) may be shared global-pY variation rather than signalling.

## 4. Alternative biology for the RTK-independent adaptor pY (`A2_adaptor_residual_correlates.tsv`)
Correlates of the residual adaptor pY (after RTK, purity and proliferation):
- **SFK substrate composite** (PTK2 Y576/Y577/Y925/Y861, PXN Y118/Y31, CAV1 Y14, BCAR1, CTTN, CDCP1): ρ = 0.36–0.60 in CCRCC, GBM, LUAD, PDAC and UCEC; 0.06–0.23 in HNSCC, LSCC and OV. The SRC/YES/FYN activation loop is not quantified.
- **Global pY:** ρ 0.1–0.6.
- **Immune/stromal (ESTIMATE):** ρ −0.25 to 0.30, mostly about 0.1. PTPRC protein: |ρ| ≤ 0.22. Purity: |ρ| ≤ 0.08.

SFK + immune + stroma + bgPY together explain 2–44% of the residual's variance (median ~0.2). Adjusting the transmission estimate for immune + stroma leaves the RTK-independent part unchanged (0.44). Adjusting for SFK or for everything gives 0.42 and 0.37. The RTK-independent pY therefore looks like an SFK / global tyrosine-phosphorylation axis, not infiltrating immune or stromal cells. That fits the literature: PTPN11 Y542/Y580 are SFK and ABL sites, and SHC is an SFK substrate.

## 5. Cell-line contrast (compendium; `A4b`, `A4`)
EGFR axis: 47 conditions from 8 studies (6 EGF / serum / HB-EGF and 41 EGFR/ERBB2 inhibitors at ≥100 nM; 33 are decryptM dose points, which are not independent).
- **Same direction of change:** RTK & adaptor 76–83% of conditions, adaptor & ERK 67–85%, RTK & ERK 80–83%.
- **Pooled correlations:** RTK→adaptor 0.77, adaptor→ERK 0.66, RTK→ERK 0.58. The product (0.51) is close to the observed RTK→ERK: transitive.
- **Decomposition:** RTK-driven part 0.79 [0.50, 1.28] vs residual 0.08. Within study: 0.66 [0.42, 0.94] vs 0.38. This is the reverse of the tumour `tc` pattern.
- **All growth-factor plus RTK-inhibitor conditions** (182 conditions, 13 studies): IV 0.43 vs residual 0.54 (difference CI −0.46 to 0.21); RTK→ERK 0.33 vs product 0.31.

**Implication.** Acute changes in receptor activity propagate. Chronic, tumour-to-tumour differences in receptor *abundance* raise receptor and adaptor pY but leave ERK output unchanged or slightly lower. That is consistent with robust steady-state MAPK output under negative feedback or adaptation (ERK→RAF/SOS feedback, DUSP/SPRY induction) and with amplified-RTK tumours routing signal elsewhere, e.g. PI3K over MAPK at high EGFRvIII. The cell-line evidence is limited: few EGF conditions, mostly decryptM, and short time scales.

## 6. Literature
- Fritsche-Guenther et al. 2011 (Mol Syst Biol 7:489): steady-state ERK phosphorylation is robust to protein-level perturbation through a single ERK→RAF1 negative feedback. This supports abundance decoupling at steady state. [link](https://doaj.org/article/89ef771c22b44fc196560c75b2a4eca5)
- Huang et al. 2007 (PNAS, EGFRvIII titration): at high EGFRvIII expression PI3K dominates over MAPK, and ERK pY correlated negatively with GBM cell growth. [PMC1937558](https://pmc.ncbi.nlm.nih.gov/articles/PMC1937558)
- Glioma proteogenomics (Nat Commun 2023, s41467-023-36005-1): EGFR amplification/mutation goes with ERK5 protein, not ERK1/2. [link](https://www.nature.com/articles/s41467-023-36005-1.pdf)
- Earlier GBM work on ERK1/2 activation vs EGFR amplification: [UV record](https://producciocientifica.uv.es/documentos/5ec07fff299952560fd24549?lang=en). I could retrieve only the title.
- SFKs phosphorylate SHP2 Y542/Y580 (Reactome R-HSA-1433488), and ABL also phosphorylates Y580. [Reactome](https://reactome.org/content/detail/R-HSA-1433488), [PMC2976378](https://pmc.ncbi.nlm.nih.gov/articles/PMC2976378), [PMC2681299](https://pmc.ncbi.nlm.nih.gov/articles/PMC2681299)
- Feedback reviews on SPRY/DUSP: [Aston review](https://publications.aston.ac.uk/id/eprint/28571/1/Negative_feedback_regulation_of_the_ERK1_2_MAPK_pathway.pdf)
- Mertins 2016 / Krug 2020 CPTAC BRCA (ERBB2 amplicon phospho): [link](https://edoc.mdc-berlin.de/id/eprint/16835/)

Novelty: "RTK amplification/overexpression does not raise ERK in tumours" is expected from these sources. New here is the quantitative, cohort-wide statement, including the negative abundance coefficient, and its absence under acute perturbation.

## Caveats of this audit
- I reused I24's site definitions so the comparison is like-for-like.
- Adaptor and ERK composites mix singly and doubly phosphorylated peptides, and in BRCA and LSCC the ERK stage is pY-only.
- RTK sites are missing not at random (low expressers drop out).
- Bootstrap SEs are IQR-based to tame weak-IV tails.
- The simulation assumes independent errors.
- Compendium conditions are not independent; decryptM dose series dominate.
