# I25 SUMMARY: is the reverse layer drug-induced redistribution of cell-cycle / culture state? Answer: **No** (all four pre-registered hypotheses unsupported)
Pre-registration: PREREG.md (written before outcomes). Test set: 285 audit_I19 inhibitor experiments from 31 studies (binders masked). CDK, CDK4_6 and AURK_PLK classes and the axis-defining studies are excluded; a sensitivity run includes those classes. Stimulation mirror: 62 experiments.
Reference axes come only from conditions outside the test set (robust z, consensus |z| ≥ 1.5):
- **M** (mitosis vs interphase): includes Olsen 2010 (20068231) and Dephoure 2008 (18669648), both already in the compendium. CC_Mup 1,405 sites; CC_Mdown 1,231.
- **MKI**: RO3306/AZD1152/MLN8054/BI2536 footprints, 4,262 sites.
- **G1/S/G2**, 4,209 sites.
- **Metformin** 24 h / 30 min (12 lines), **AMPK/MARK/SIK motif**, **3D vs 2D**, **Hippo genes**, **serum**, **OIS**.

Controls: strata of residue × coverage × variance (mitotic sites have 1.8× higher compendium SD); matched random removal (20 draws); study-cluster bootstrap.

| Hypothesis | Key estimate [95% CI] | Verdict |
|---|---|---|
| **H1**: reverse responders enriched for CC sites | CC_Mup log(O/E): reverse **−0.12 [−0.32, 0.08]**, primary −0.33 [−0.59, −0.09]. CC_union: reverse +0.15 [0.05, 0.25] vs primary +0.08 [−0.01, 0.16]; reverse−primary **0.06 [−0.06, 0.20]**. In Hijazi (165 exps) O/E ≈ 1.05 in both directions. Unmatched results are the same (T1b). With CDK/AURK classes included: 0.03 [−0.09, 0.16]. | Not supported |
| **H2**: removing CC sites restores reverse replication | Cross-study same class (\|log2FC\| > 1): LOR_up 0.40 vs LOR_down 0.83, gap 0.43. Minus CC_union: up 0.39 (Δ −0.005; matched random Δ +0.015); **gap closure ≈ 0%** (Mup/MKI/Mdown/inter: −4% to +6%). Same cell line across studies (incl. MCF7 Hijazi vs Ochoa): up 0.61 → 0.66 (+0.05; random +0.01, 97.5th percentile +0.07); closure 9%. Within CC sites, reverse changes replicate *better* (LOR_up 0.65) than reverse changes overall (0.40). | Not supported; CC sites are not the discordant ones |
| **H3**: within study, reverse magnitude tracks the cell-cycle shift; the shift is lab-specific | Partial ρ(shift_M, reverse magnitude on non-CC sites \| primary magnitude, frac_down) = **0.04 [−0.21, 0.30]**; 6 units, 223 exps, 3/6 positive. Shift_M is **portable**: same drug in a different cell line (Hijazi) r 0.57; same class across studies r 0.57 (12 unit pairs); MCF7 Hijazi vs Ochoa same drug r 0.32 (n 11). Variance (units with ≥3 exps): shift_M is **class 0.41 / unit 0.07**, whereas reverse magnitude is **unit 0.50 / class 0.09**. Caveat: primary magnitude is also unit 0.50. | Not supported, and inverted: the cell-cycle shift behaves like the primary layer |
| **H4**: other culture axes or normalisation | Site removal: no axis closes ≥ 50% of the gap. Best is metformin-24 h at 12%, and that comes from LOR_down falling, not LOR_up rising; Hippo, serum, OIS and 3D are ≤ 4%. Enrichment: none for reverse-specific sites. Hippo sites are depleted among *primary* responders (−0.43 [−0.68, −0.20]). Exp-level: only AMPK/MARK/SIK motif shift, partial ρ 0.22 [0.03, 0.42] (nominal, 1 of 12 features). Pair-level variance in reverse replication (1,081 cross-study pairs): every axis has R² < 0.03; \|Δfrac_down\| 0.08 [0.00, 0.28], **with the opposite sign** (r +0.29); joint 12 predictors 0.085; same cell line 0.003; \|Δtime\| 0.02. Per-experiment median centring does not close the gap (0.43 → 0.48). | Not supported; normalisation excluded |

Interpretation:
- The reverse layer's lab/culture specificity is **not** explained by an externally defined cell-cycle, energy, density/Hippo, serum, senescence or normalisation axis.
- Drugs do shift the mitotic axis, but in a reproducible, class-specific way (part of portable pharmacology).
- What remains is a study×culture-specific component that sits orthogonal to all the state axes tested.
- Prime remaining candidate (not tested here): shared-control / within-run structure. I22's centring control addressed this only in Hijazi.
- Incidental findings in the stimulation mirror:
  - Stimulus primary (up) responders are enriched for metformin-24h-up sites (+0.26 [0.15, 0.37]) and the stimulus reverse layer for metformin-down sites (+0.29 [0.11, 0.47]). This is a shared stress-signalling footprint, not reverse-specific.
  - The stimulus reverse (down) layer is mildly enriched for CC_union: 0.25 vs primary 0.07, diff 0.18 [0.02, 0.34]. This is the only CC-positive result. It is uncorrected for the 14 sets tested, so treat it as a lead at most.

Novelty:
- Olsen 2010 showed that most phosphosites are cell-cycle regulated.
- Liu et al. (HeLa heterogeneity across 13 labs, bioRxiv 307421 / Nat Biotech 2019) shows cross-lab biological heterogeneity, but not in phosphosite drug responses.
- I found no explicit test of cell-cycle confounding of inhibitor compensatory responses (15 min search). This is a modest negative result that prunes a plausible mechanism.

**Self-score: Compelling 3** (multiple concordant, controlled negatives; limited by HeLa-centric, partly low-coverage reference axes) **/ Novel 2 / Promising 2.**

Negatives and caveats:
- The axes are HeLa-heavy; the Olsen ratios are on an odd scale (robust z handles this); only 1,085 CC_Mup sites fall in the test universe.
- Most treatments are ≤ 1 h, so population redistribution is a priori unlikely. The time strata (> 60 min, 4 studies) are underpowered.
- R² bootstrap lower bounds are ≥ 0 by construction, so the R² ≥ 0.10 threshold is the decision rule.
- Deviations from PREREG: +0.5 pseudocount (zero counts); added the unmatched control (T1b), the units-≥3 variance table (T3d) and the normalisation test (T2b).

Files:
- `scripts/b_I25/01_axes.py`, `common_i25.py`, `02_h1.py`, `02b_h1_agg.py`, `03_h2.py`, `03b_h4_norm.py`, `04_h3.py`, `05_h4_pairs.py`, `06_figure.py`
- `results/b_I25/tables/`: T1, T1b, T2, T2b, T3, T3b, T3c, T3d, T4
- `figures/fig1_reverse_layer_state_axes.png`
- Axis cache: `cache/b_I25/axes.parquet`
