# Audit of I19 key finding: up/down asymmetry in cross-study replication under kinase inhibition

Auditor scripts: `/workspace/scripts/audit_I19/` (a01–a08). Outputs: `/workspace/results/audit_I19/` (T_*.tsv, fig_mirror_inhibitor_vs_stimulation.png). b_I19 was not modified.
Primary audit metric (my own implementation): experiment pairs (study × cell × drug) from **different studies**, same disjoint inhibitor class (23 classes with ≥2 studies, 41 studies); raw log2FC vs control; own-binder proteins masked; sites measured in both. Discovery = dir·x_A in |log2FC| bin. Replication = dir·x_B > 0.3. Reported as **log odds ratio (LOR)** of replication for discovery vs non-discovery co-measured sites in B. This is the hypergeometric-type null and matches per-direction, per-experiment base rates. Each study pair is weighted equally, and CIs come from a cluster bootstrap over studies.

## Verdict
The claim **survives with qualification.** The core asymmetry is real and is not a measurement artefact. The best evidence is the mirror test: in activating stimulations, increases replicate and decreases do not. Several parts of the wording are overstated, though:
- "Increases do not replicate" is false. They replicate weakly but above the base rate (LOR ≈ 0.3, OR ≈ 1.4).
- The inhibitor asymmetry is modest and depends on the MEK class.
- Most of the "portable down" signal is not specific to the drug class.
- The decomposition does not show that the asymmetry is about the study changing. Both directions lose replication when the cell line changes; decreases lose less.

## 1. Code review and recompute
- **Reproduced exactly** with an independent implementation (a06) of the I19 metric: down 0.204/0.292/0.394 and up 0.125/0.134/0.143.
- **Bug B1 (null and observed averaged over different pair sets).** In `08_replication.py`, a null exists only when study B has at least as many other-group experiments as group experiments. That holds for 54% of pairs. `rep` is averaged over all pairs, but `rep_null` only over that subset. Restricted to the matched pairs (20 null draws instead of 1):

  | | Down | Up |
  |---|---|---|
  | rep | 0.168 / 0.257 / 0.393 | 0.097 / 0.113 / 0.129 |
  | null | 0.126 / 0.174 / 0.239 | 0.086 / 0.097 / 0.108 |
  | excess | +0.04 / +0.08 / +0.15 | +0.01 / +0.02 / +0.02 |

  The asymmetry survives, but the lowest-bin down excess halves.
- **B2 (double counting).** Groups overlap: AXIS_RTK, AXIS_MAPK and AXIS_PI3K contain MEK, EGFR, MET_ALK and others again. AXIS_RTK alone supplies 34% of discovery sites. The effect on the results is minor, since disjoint groups give 0.200/0.285/0.388 vs 0.129/0.142/0.154.
- **B3 (null is a single random draw).** The null uses one draw (seed 3) per study. This is noisy but unbiased.
- **Conditioning on measurement.** The metric does condition on the site being measured in both studies.
- **Centring.** It does not drive the result: raw and centred results agree within 0.02.
- **Experiment-level, base-rate-matched recompute** (a02/a03; T_headline_recompute.tsv, T_tests_a_d.tsv; bins 0.5–1 / 1–2 / >2):

  | Direction | rep | base | site-matched other-class null | LOR |
  |---|---|---|---|---|
  | Down | 0.27 / 0.32 / 0.35 | 0.22 | 0.24 / 0.28 / 0.30 | 0.32 / 0.63 / 0.77 |
  | Up | 0.20 / 0.21 / 0.20 | 0.17 | 0.19 / 0.19 / 0.19 | 0.23 / 0.34 / 0.35 |

  dLOR(down − up) is 0.30 [0.05, 0.69] for the 1–2 bin and 0.27 [−0.08, 0.87] for >2.
- **Up sites are variable sites rather than replicating sites.** When a site rises in A, its mean in B is about 0 (−0.02 to +0.01). It is also slightly more likely than the base rate to fall in B (0.22–0.24 vs 0.21–0.22).
- **Class specificity is small.** Replication above the other-class null (rep − snull) is +0.03/+0.05/+0.08 for decreases and 0.00/+0.02/+0.01 for increases. The difference has CIs that include 0 (e.g. 1–2 bin: 0.036 [−0.02, 0.12]). Most cross-study down replication is shared with other inhibitor classes in study B, i.e. convergent downstream nodes.

## 2. Alternative explanations

**a. Missingness / detection.**
- Missing in B given a discovery site in A: down 0.70, up 0.72, all sites 0.66. There is no MNAR advantage for decreases.
- Restricting to sites measured in ≥50% of the 310 experiments: down LOR 0.66/0.88 vs up 0.43/0.34 (1–2 and >2 bins).
- Restricting to the top quartile of compendium coverage: 0.65/0.86 vs 0.36/0.40, dLOR 0.29 [0.03, 0.71].
- Restricting to sites measured in all experiments of the class: underpowered. The 1–2 bin gives 0.88 vs 0.46; the >2 bin reverses (0.61 vs 0.69).
- Imputation signature: one study (19651622) has 16% of values at −0.737, consistent with imputation. It contributes very little.
- **Conclusion: not explained by missingness.**

**b. Base rates.**
- Down base rates are higher (0.22 vs 0.17). The LOR null accounts for this, and the asymmetry persists (above).
- The I19 other-class null is site-matched and is valid apart from bug B1.

**c. Time and dose matching** (condition level, T_tests_c_e.tsv). Values are down vs up LOR for the 1–2 / >2 bins.

| Matching | Down LOR | Up LOR | dLOR | Note |
|---|---|---|---|---|
| Both ≤30 min | 0.31 / 0.61 | 0.41 / 0.78 | ≈ 0 | **no asymmetry**; 10 study pairs |
| Both 31–119 min (mostly 60 min) | 0.83 / 1.16 | 0.36 / 0.18 | 0.48 [−0.44, 1.74] | 14 study pairs |
| Both ≥2 h | 0.22 / 0.40 | 0.36 / 0.30 | ≈ 0 | 5 study pairs |
| Both dose ≥1 µM | 0.57 / 0.73 | 0.34 / 0.37 | 0.24 [0.04, 0.54] | |
| Both ≤60 min and ≥1 µM | 0.54 / 0.81 | 0.38 / 0.31 | 0.18 [−0.32, 0.92] | |

- Time-matched data are thin.
- The asymmetry is clear only in the 60-min stratum (Hijazi, decryptM, Ochoa 26060313), and it is absent at ≤30 min.
- **Time mismatch does not create the asymmetry, but the data cannot show that it holds at every time point.**

**d. Study composition.**

| Subset | dLOR 1–2 | dLOR >2 |
|---|---|---|
| Drop decryptM | 0.26 [0.02, 0.63] | 0.26 |
| Drop Hijazi | 0.43 [0.15, 0.73] | 0.47 [0.03, 0.91] |
| Classes with ≥4 studies | 0.30 [0.05, 0.65] | 0.25 |
| Leave-one-study-out (range) | 0.25 – 0.43 | |

- **Leave-one-class-out: dropping MEK roughly halves the effect** (0.20 [−0.05, 0.53] / 0.13 [−0.17, 0.64]).
- Per class, the direction is positive in 9 of 12 classes with ≥4 studies, but significant only for MEK (0.70 [0.34, 1.17]).
- **DDR inhibitors show the reverse at >2:** increases replicate, LOR 1.03 vs 0.10. Likely cause: ATR/ATM inhibition leads to compensatory activation of DNA-PK/ATM, a portable compensatory response.
- SRC_ABL and mTORcat trend reversed.

**e. Noise ceiling and decomposition** (within study, condition level; LOR down / up for the 1–2 and >2 bins):

| Pair type | Down | Up |
|---|---|---|
| Same drug, same cell, different dose/time | 1.99 / 2.15 | 2.44 / 2.76 |
| Same class, different drug, same cell | 2.25 / 2.89 | 1.90 / 2.27 |
| NULL: different class, same cell (shared control/batch) | 1.24 / 1.64 | 1.22 / 1.27 |
| Same drug, different cell (same study) | 1.25 / 1.35 | 0.80 / 0.92 |
| Cross-study, same cell | 1.06 / 1.20 | 0.58 / 0.43 |
| Cross-study, different cell | 0.53 / 0.67 | 0.31 / 0.36 |

- **Increases are not noisier.** Within the same drug and cell, they replicate as well as or better than decreases (dLOR −0.45/−0.75).
- **Changing the cell line is the step that opens the gap.** Same drug, different cell gives dLOR 0.45 [0.19, 0.96], even within one study.
- Both directions lose replication at each step away from the same context; decreases lose less.
- I19's statement that changing cell line or study "destroys" replication for increases is directionally right. But I19's own pair-type table shows decreases also fall to near-null on its median-excess metric (+0.09 / +0.03).

**f. Mirror test in non-inhibitor perturbations** (key discriminating test; a05, T_tests_f_stim.tsv)
- Data: 75 experiments, 52 studies. Pooled activating classes are GF/RTK (13 studies), TNF (6), DNA damage (12), GPCR (7), phosphatase inhibitors (4) and PMA. Same metric as for inhibitors.
- Pooled activating classes:
  - Increases: LOR 0.56 / 1.11 / 1.62 (rep 0.30 / 0.40 / 0.50 vs base 0.20).
  - Decreases: LOR 0.48 / 0.66 / 0.59 (flat).
  - dLOR(down − up) is −0.38 [−0.68, −0.16] for 1–2 and −0.70 [−1.21, −0.25] for >2.
- By class (sign of dLOR at 1–2 / >2):

  | Class | 1–2 | >2 |
  |---|---|---|
  | GF/RTK | −0.57 | −0.48 |
  | GF/RTK, ≤60 min only | −0.63 | −0.58 |
  | TNF | −1.06 | −0.89 |
  | DNA damage | −0.22 | −0.89 |
  | GPCR | −0.09 | −0.40 |
  | Phosphatase inhibitors | +0.13 (n.s., 4 studies) | +0.27 (n.s.) |

- Robust to leave-one-class-out (−0.26 to −0.49) and leave-one-study-out (−0.33 to −0.43).
- **The mirror image is present and stronger than in inhibitors.** In both settings the primary direction of the perturbation is portable and the reverse direction is not. This rules out a generic technical bias favouring decreases (missingness, quantification floor, base rates) and supports "primary-direction effects are portable; reverse/secondary effects are context-specific".

## 3. Defensible wording
> "Across 41 studies, phosphosite changes in the perturbation's primary direction replicate across studies more strongly, and more steeply with effect size, than changes in the opposite direction.
> - For kinase inhibitors, decreases have cross-study replication LOR 0.32 → 0.63 → 0.77 over |log2FC| bins, vs 0.23 → 0.34 → 0.35 for increases (Δ 0.30 [0.05, 0.69] at |log2FC| 1–2). The effect is strongest for MEK inhibitors and attenuated without them.
> - Activating stimulations show the mirror image (increases 0.56 → 1.11 → 1.62 vs decreases 0.48 → 0.66 → 0.59).
>
> Reverse-direction changes are reproducible within the same cell line (as reproducible as primary changes) but lose more replication when the cell line changes. This is consistent with secondary, compensatory responses being more cell-context-dependent. Most cross-study-portable decreases are shared across inhibitor classes, i.e. convergent pathway nodes rather than drug-specific footprints."

Avoid these phrasings:
- "increases do not replicate" (they do, weakly)
- "direct inhibition footprints" (not shown to be direct or class-specific)
- "lost when the study changes" (the cell line, not the study, is the main step)
- quoting the 0.12/0.16/0.22 null (bug B1)

## 4. Most compelling single figure
`fig_mirror_inhibitor_vs_stimulation.png` (data: T_fig_mirror.tsv) shows cross-study replication LOR vs |log2FC| bin for increases and decreases, with kinase inhibitors and activating stimulations side by side and study-cluster bootstrap 95% CIs. In inhibitors the decrease curve rises above the flat increase curve; in stimulations the curves swap. A within-cell vs cross-cell panel (table e) would be the natural second panel.

## Self-assessment of the claim after audit
- Confidence: moderate (C3). The inhibitor-only CIs are wide and MEK-dependent; the mirror test is robust.
- Novelty: N3–N4. The generalised "primary-direction portability" framing is more novel than the inhibitor-only claim.
