# val_B: does measured cold-ischaemia time explain the candidate-B MAPK-axis coherence?

Scripts: `/workspace/scripts/val_B/01_fetch.py` (GDC), `02_table.py`, `03_analysis.py`, `04_meta_erk.py`. Tables: `/workspace/results/val_B/` (`T1_coverage`, `T2_correlations`, `T3_coherence`, `T4_meta_time_corr`, `T5_erk_vs_isch_adjust`, `gdc_ischaemia_all.tsv`). Raw GDC JSON: `/workspace/cache/val_B/dl/`.

## Bottom line
In the 7 cohorts that have measured times (n = 723 tumours), the footprint ischaemia score does not track measured excision-to-freezing time. Regressing measured time (linear + log) out of every site leaves MEK/EGF/TNF/DNA-damage coherence unchanged (within about 2%, same as a permuted-time control). Regressing the footprint ischaemia score still halves MEK/EGF/TNF coherence. So the "ischaemia footprint" is a general MAPK/stress-output axis readout, not a measure of pre-analytical time, and the audit's "may be pre-analytical" inference is not supported by measured data. Measured time is narrow-range and covers 7/10 cohorts, so some other pre-analytical effect cannot be excluded.

## 1. Coverage (T1; GDC `samples.time_between_excision_and_freezing`, tumour samples, median per case; CPTAC-3 and CPTAC-2 projects queried)
| cohort | tumours (phospho) | with excision time | median [IQR] min | range | clamping time (n) |
|---|---|---|---|---|---|
| CCRCC | 103 | 103 (100%) | 17 [12-23] | 2-30 | 25 |
| GBM | 99 | 99 (100%) | 12 [8.5-15] | 2-29 | 0 |
| HNSCC | 108 | 108 (100%) | 10 [8-14] | 3-30 | 14 |
| LSCC | 108 | 107 (99%) | 15 [10-22] | 1-425 | 56 |
| LUAD | 110 | 106 (96%) | 13 [7-21] | 1-40 | 70 |
| PDAC | 105 | 105 (100%) | 20 [10-35] | 5-50 | 21 |
| UCEC | 95 | 95 (100%) | 16 [12-22] | 3-30 | 1 |
| **BRCA, COAD, OV** | 122 / 97 / 83 | **0** | n/a | n/a | 0 |

- BRCA, COAD and OV are CPTAC-2 (TCGA-era) cases. Their GDC records have the field empty. A scan of the GDC case field list found no other ischaemia/freezing/warm-time field. These 3 cohorts remain unmeasured.
- One LSCC value (525,618 min) is an obvious entry error. I set it to NA (cut-off >1440 min). LSCC then has 107 cases.
- The measured window is short (IQR roughly 8-35 min; max 50 except one LSCC case at 425 min). It reflects protocol-compliant collection. The power is limited accordingly (see caveats).
- Clamping time (warm plus cold) is too sparse to use.

## 2. Correlations with measured time (Spearman, per cohort; T2; pooled in T4)
| variable | CCRCC | GBM | HNSCC | LSCC | LUAD | PDAC | UCEC | pooled rho [95% CI] (fixed) |
|---|---|---|---|---|---|---|---|---|
| footprint ischaemia score | .06 | .31 | -.01 | -.05 | -.16 | -.20 | -.16 | **-0.03 [-0.11, 0.04]** |
| MEK score | .15 | -.16 | .18 | .09 | -.05 | .30 | .08 | 0.09 [0.01, 0.16] |
| EGF score | -.04 | .20 | -.09 | -.01 | .17 | -.24 | -.11 | -0.02 [-0.09, 0.06] |
| TNF score | .09 | .19 | -.07 | -.04 | .01 | -.11 | -.23 | -0.02 [-0.10, 0.05] |
| DNA-damage score | .07 | .20 | .14 | .04 | -.09 | .31 | .31 | 0.14 [0.07, 0.21] |
| ERK pTpY (T185/Y187/T202/Y204 mean z) | -.10 | .20 | -.17 | -.07 | .11 | -.39 | .10 | -0.05 [-0.13, 0.02] |

- The LUAD result is reproduced (ischaemia score rho -.16). It is not specific to LUAD: signs differ across cohorts, with heterogeneity I2 about 0.7, and the pooled estimate is near 0. The CI excludes rho of about 0.11 or more in either direction.
- Time explains 0.3-8% of the variance of the footprint ischaemia score (R2 on t + log t, per cohort).
- If the score were an ischaemia readout, longer time should raise EGF/TNF scores and lower the MEK score, as the score's own correlations imply (EGF +.4 to +.8, MEK -.6 to -.8 in each cohort). The observed pattern is not that: MEK is weakly positive, and EGF/TNF are about 0.
- DNA damage is the only footprint with a (weak) positive pooled rho. That is the footprint the audit said is unaffected by the ischaemia score, and it is heterogeneous (-.09 to .31). I do not interpret it.
- ERK pTpY does not rise with ischaemia time either (pooled -0.05; PDAC -.39 is the one notable cohort, and it is negative).
- The footprint ischaemia score itself correlates with MEK (-.56 to -.76), EGF (+.38 to +.81), TNF (+.47 to +.80) and ERK pTpY (+.15 to +.44) in every cohort. These values are unchanged after partialling out log time (partial rho within .03 of the raw). So the score tracks the MAPK/stress axis but not time.

## 3. Coherence after adjustment (T3; median signed coherence over cohorts with times: 6-7 cohorts; each adjustment is applied to the same time-known tumours)
| footprint | none | measured t + log t | measured log t | permuted t (control) | footprint ischaemia score | footprint score + t + log t | val_E isch (pathway sites removed) |
|---|---|---|---|---|---|---|---|
| MEK | .170 | **.169** | .167 | .169-.173 | **.087** (-49%) | .086 | .088 |
| EGF | .145 | **.145** | .145 | .144-.146 | **.082** (-43%) | .080 | .082 |
| TNF | .073 | **.072** | .071 | .072-.073 | **.028** (-62%) | .031 | .027 |
| DNA damage | .071 | .068 | .067 | .066-.071 | .066 (-7%) | .064 | .066 |

- Measured time does not reduce coherence. In no cohort-set pair does it change coherence by more than 8% (largest fall UCEC TNF -8%; largest rise +4%). The effect is indistinguishable from permuted time.
- The footprint ischaemia score reproduces the audit (all-tumour, 10 cohorts: MEK .153 to .084, EGF .141 to .074, TNF .066 to .025, DNA damage .064 to .065). Adding measured time on top does nothing further.
- Removing the pathway-protein sites from the score (val_E `isch_noPath`) does not change its effect. The score's reach is therefore not due to a handful of MAPK sites. It sits on the broad phospho axis shared by the MEK/EGF/TNF footprints (the audit's top-PC axis).

## 4. Is the ischaemia footprint "another MAPK readout"? (T5, all 10 cohorts)
Regressing an independent ERK activation-loop readout (T185/Y187/Y204 mean z; T202 left out because it is in the ischaemia set) from every site gives median coherence MEK .153 to .102 (-33%), EGF .141 to .114 (-19%), TNF .066 to .059, DNA damage .064 to .065. The ischaemia score does more than that (.084 / .074 / .025). After conditioning on the ERK readout, the ischaemia score still correlates with MEK, EGF and TNF scores (median partial rho -0.62, +0.64, +0.75). It is therefore a broad "stimulus/stress-responsive phospho" axis (the audit's PC1), not only ERK activity. It overlaps MAPK output (rho with ERK readout .34) and is not a measure of cold-ischaemia time.

## Verdict
1. **Measured ischaemia does not explain the MAPK-axis coherence.** Not with the correlations (pooled rho of time with MEK/EGF/TNF scores is between -0.02 and +0.09) and not with the adjustment (coherence unchanged). The audit's halving came from adjusting for a footprint score that is a readout of the shared MAPK/stress axis. It does not track measured time in 7 cohorts (pooled rho -0.03, no cohort-consistent sign).
2. The compendium "ischaemia footprint" (Mertins 2014 top sites) is therefore not a validated pre-analytical proxy. Using it as a covariate over-adjusts, because it removes real axis variance. Audit finding B3 should be downgraded from "major confounder" to "unresolved shared axis, not explained by measured cold-ischaemia time". This also affects the pervanadate and GPCR statements that relied on it.
3. **Caveats:**
   - 3/10 cohorts (BRCA, COAD, OV) have no measured time.
   - The measured range is narrow (about 8-35 min IQR). A confound acting only above about 1 h would not be seen.
   - Warm ischaemia, time to resection, anaesthesia and tissue handling are not captured. The pooled CI (about +-0.07) rules out time effects larger than about rho 0.1, not smaller ones.
   - The same axis could still be driven by an unmeasured pre-analytical factor, but it is not driven by recorded cold-ischaemia time.

## Defensible wording for candidate B
"Consensus perturbation footprints from cell-line phosphoproteomics remain coherent in all 10 CPTAC cohorts, with the MEKi/EGF/TNF footprints forming a single MAPK/stress-output axis that tracks ERK activity. This axis does not track recorded cold-ischaemia time: in the 7 cohorts with GDC excision-to-freezing times (n = 723, median 10-20 min), time is not correlated with the footprint scores (pooled |rho| <= 0.09) and regressing it out leaves coherence unchanged. A compendium-derived ischaemia signature correlates strongly with the axis, but it does not track measured time either, so it is a co-readout of the axis and not evidence of a pre-analytical origin. Collection effects outside the recorded cold-ischaemia window (3 cohorts without times, warm ischaemia, handling) cannot be excluded."
