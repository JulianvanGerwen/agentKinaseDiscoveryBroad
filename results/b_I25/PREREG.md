# I25 PREREG — is the reverse layer drug-induced redistribution of cell-cycle / culture state?
Written 2026-10-29 BEFORE any outcome computation. Only condition metadata and reference-contrast orientation/scale were inspected.

## Data
- Test set: audit_I19 inhibitor experiments (`results/audit_I19/cache/inh_exp.parquet`, study×cell×drug, sign-corrected, own binders masked).
  Primary direction = down; reverse = up. **Excluded from test set**: classes CDK, CDK4_6, AURK_PLK (their primary effect is on cell-cycle sites)
  and all studies used to define reference axes (21712546, 21857030, 26055452, 26270265, 27134283, 29970458). Sensitivity: including CDK/AURK classes.
- Mirror (secondary): audit stimulation experiments (reverse = down).

## Reference axes (defined only from conditions outside the test set)
Per reference condition: robust z (x − median)/(1.4826·MAD). Site score = mean z over the axis contrasts in which it is measured.
- **M axis (mitosis vs interphase; includes Olsen 2010 = PMID 20068231 and Dephoure 2008 = 18669648, both already in the compendium):**
  20068231_qPTM_2 (mitosis; oriented M-up as judged from CDK1 Y15↓/NEK9 S869↑), 18669648_Ochoa_267 (nocodazole M vs async), 18691976_qPTM_1 (M vs S),
  28674419_qPTM_2 (mitotic vs control), 29052541_qPTM_2 (M vs G1, NB4).
- **MKI axis (mitotic-kinase inhibitor footprint, sign flipped so + = mitotic-kinase dependent):** RO3306 (26055452_1, 27134283_2), AZD1152 (26270265_1),
  21712546 MLN8054 ×2, AZD1152+ZM, BI2536 ×2.
- **Interphase axis (G1/S/G2; secondary):** 32051232_1–3, 29052541_1/_3, 23322592_1, 18669648_Ochoa_72, 29154753_1 (|z| based).
- Sets: **CC_Mup** = M-score ≥ 1.5 and ≥75% of contributing contrasts z>0 (primary set); CC_Mdown = M-score ≤ −1.5 (same consistency);
  CC_MKI = MKI-score ≥ 1.5; **CC_union** = CC_Mup ∪ CC_Mdown ∪ CC_MKI (used for H2 removal); CC_interphase = |interphase score| ≥ 1.5.
- Universe for enrichment = sites measured in ≥1 contrast of the axis (avoids reference-coverage confounding).

## Competing axes (H4)
- ENERGY: metformin 24 h and 30 min across 12 CRC lines (36781298): consensus |mean z| ≥ 1.5, ≥75% sign agreement (24 h primary). Also motif set AMPK/MARK/SIK/NUAK/BRSK (Cantley percentile ≥ 95 for any of PRKAA1/2, MARK1–4, SIK1–3, NUAK1/2, BRSK1/2).
- CULTURE/DENSITY proxy: 3D vs 2D (30061712, |mean z| ≥ 1.5); HIPPO gene set (all sites on YAP1, WWTR1, LATS1/2, STK3/4, MOB1A/B, SAV1, NF2, AMOT, AMOTL1/2, WWC1, PTPN14, FRMD6, TEAD1–4).
- GROWTH/SERUM: serum time course 2–48 h (29970458) and OIS senescence (24961811) consensus |mean z| ≥ 1.5.
- NORMALISATION: experiment-level fraction of sites with x < −0.3 (frac_down) and median x; used in H3 and pair models.

## Controls
- Matched background: strata residue (S/T/Y) × compendium coverage decile (n_conditions) × variance decile (SD across non-test conditions); expected set membership of responders computed from stratum-specific rates within each experiment.
- Matched random removal (H2): 20 draws of non-axis sites with the same stratum composition as the removed set.
- Study-cluster bootstrap (1000 for H1, 300–500 for pair metrics) for all CIs; LOSO for main effects.

## Hypotheses, statistics, and decision rules
- **H1**: log(O/E) of CC_Mup membership among reverse responders (x > 1) vs matched background, aggregated per study then averaged; compare with primary responders (x < −1).
  Supported if reverse log(O/E) > 0 with 95% CI excluding 0 AND reverse − primary > 0 with CI excluding 0. Also reported for CC_union, CC_Mdown, CC_MKI, CC_interphase.
- **H2**: audit_I19 cross-study pair metric (same disjoint class, different study, LOR of replication, bins 0.5–1/1–2/>2; summary = 1–2 bin and pooled >1), computed on (a) all reference-covered sites, (b) excluding CC_union, (c) excluding matched random sets; also within CC_union only.
  Supported if removal raises LOR_up by ≥ 50% of the (down − up) gap AND the rise exceeds the matched random-removal 95% range. Secondary: same-cell cross-study pairs (incl. MCF7 Hijazi vs Ochoa).
- **H3**: per experiment, shift_M = mean x over CC_Mup − mean x over CC_Mdown (needs ≥10 sites each); reverse magnitude R = fraction of non-CC_union sites with x > 1.
  Within study×cell units with ≥6 experiments: Spearman ρ(shift_M, R), and partial ρ controlling for primary magnitude (fraction x < −1) and frac_down; averaged across units, unit bootstrap.
  Second part: for the same drug across cell lines (Hijazi) and same class across studies, concordance of shift_M (Pearson r) vs concordance of primary footprint magnitude; variance decomposition of shift_M into unit vs class.
  Supported if mean partial ρ > 0.2 with CI excluding 0 AND cross-unit concordance of shift_M < within-unit (same drug different dose/time) and < 0.3.
- **H4**: same three frameworks for each competing axis. Variance in cross-lab discordance: per cross-study same-class pair, outcome = LOR_up (pooled >1 bin); predictors = |Δ axis shift| between the two experiments (for M, MKI, ENERGY, CULTURE, SERUM, HIPPO) and |Δ frac_down|; report univariate and joint R² with pair-cluster (study) bootstrap. Also site-level gap closure per axis.
  An axis "explains" discordance if gap closure ≥ 50% or R² ≥ 0.10 with CI excluding 0.

## Expected failure modes (declared)
Reference axes are mostly HeLa and partially low coverage; mitotic sites are high-variance sites (hence variance matching); 1-h drug treatments cannot change cell-cycle distribution much (strong prior against H1 for ≤1 h; a positive H1 would then point to kinase-level coupling rather than population redistribution — interpret accordingly by stratifying by time ≤60 min vs >60 min).
