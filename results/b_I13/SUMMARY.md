# I13: Do phosphatase lesions in tumours shift sites that match the phosphatase's sequence rules? (2026-10-29)

**Verdict: mostly negative.** No phosphatase *genetic* lesion (PPP2R2A CNV loss, PPP2R1A hotspot mutation, PPM1D gain, loss of a PTP) shifts the sites that match that phosphatase's sequence rules beyond what random-gene nulls produce. The one replicated association is with low PPP2R2A *protein*, and it is not genetically anchored. The pipeline does detect a known phosphatase–substrate effect (PTEN loss raising PRAS40 pT246), so the null is informative rather than a sign of a broken test. Scripts: `scripts/b_I13/01–07`. Tables: `tables/`. Figures: `figures/fig1–4`.

**Literature (novelty check).** No published study tests phosphatase sequence rules against phosphatase genetic lesions in tumour phosphoproteomes. Dou 2020 (CPTAC UCEC) reports no PPP2R1A-specific phosphosite analysis. The PPP2R2A-loss phosphoproteome is known only from a BT474 knockdown (HUPO 2019 abstract, 56 sites up, no motif analysis). Haesen 2016 showed that PPP2R1A mutants are dominant-negative and raise GSK3β/AKT/S6K phosphorylation in cell lines. Kauko 2020 found that PP2A inhibition preferentially affects pThr. Gräf 2022 found PPM1D targets are SQ/TQ sites with flanking Glu. The rules themselves come from Cundell 2016 and Kruse 2020: B55 recognises basic residues on both sides of S/T-P sites and prefers pT.

## 1. Perturbation side (compendium): not enough data to learn rules
- Phosphatase perturbations in the compendium:
  - PP2A inhibitors: okadaic acid (OA) and cantharidic acid (CA), MCF7, 1 h, one study (26060313).
  - Pervanadate: 3 conditions from 2 studies.
  - SHP099 dose series; PRL3i.
  - Nothing for calyculin, calcineurin inhibitors or phosphatase knockdowns. External downloads (Kruse 2020, Kauko 2020) were blocked by the network.
- **OA and CA disagree** (r = 0.21 across sites). Held-out AUC (train on one, test on the other) is 0.45–0.49 for every feature set, so the model does not transfer.
- CA alone fits B55 expectations. Per-SD log-odds of a site being up (z ≥ 2): basic residues on both sides +0.39 (p = 1e-24), pT +0.10 (p = 0.02), acidic −0.30 (p = 1e-9). CA ranks at the 98.6th percentile of 1,798 compendium conditions for enrichment of bipartite-basic sites.
- OA shows the opposite pattern: +1 Pro −0.34, +1 hydrophobic +0.40. Its bipartite-basic enrichment sits at the 11th percentile.
- **Pervanadate**, leave-one-study-out AUC: nuisance 0.65–0.70; kinase-motif (multi-hot of top-3 predicted kinases) 0.70–0.77; phosphatase features 0.70–0.77; both combined 0.70–0.78. Phosphatase features add ≤0.01 beyond kinase motifs, and pY drives almost all of it. Protein-level SLiM distances (LxxIxE, RVxF, PxIxIT, LxVP) had no consistent effect.

## 2. Tumour natural experiments (CPTAC, 10 cohorts, 1,030 tumours)
**Method.**
- Stage 1, per cohort: protein-adjusted site, z-scored, regressed on the exposure, adjusting for purity, E2F/G2M Hallmark, CNV burden, ESTIMATE immune/stroma, TP53 mutation and 5 proteome PCs.
- Stage 2: the per-site slopes regressed on sequence features, with site SD, mean and observed fraction as covariates.
- Inference uses empirical nulls built from the same pipeline: 300 random CNV genes on other chromosomes, about 60 genes on 8p, 1,000 random proteins (the same proteins in every cohort), and 300 label permutations for the hotspot. Cohorts are combined by DerSimonian–Laird meta-analysis using the null SD as the SE.

**First stage: PPP2R2A CNV → protein holds.** Spearman ρ = 0.24–0.69 across cohorts. Shallow deletion lowers PPP2R2A protein by 0.17–0.55 log2.

| exposure → feature (pS/pT sites) | RE β [95% CI] | cohorts positive | verdict |
|---|---|---|---|
| (a) PPP2R2A CNV loss → B55 literature score | 0.0001 [−0.003, 0.004] | 6/10 | null |
| PPP2R2A CNV loss → bipartite-basic sites | 0.0024 [−0.0005, 0.005] | 8/10 | ns |
| PPP2R2A CNV loss → pT | 0.0004 [−0.001, 0.002] | 4/10 | null |
| PPP2R2A protein low → bipartite-basic sites | **0.0068 [0.0036, 0.0099]**, FDR 0.008 | **10/10** | replicates, but see below |
| PPP2R2A protein low → +1 Pro sites | −0.006 [−0.011, −0.001] | 2/10 | opposite to the B55 rule |
| (b) PPP2R1A hotspot (UCEC only, n = 8) → B56 score, B55 score or pT | permutation z 0.3–1.5 | — | underpowered null |
| (c) PPM1D gain → SQ/TQ sites, SQ+acidic sites; H2AX S139 | z ≈ 1 (wrong direction) | 5–6/10 | null. p53 S15 and CHK1 S345 are below the coverage threshold |
| (d) PTPRD, PTPN1, PTPN11 or PTPN12 CNV loss → pY | all \|z\| < 1 | 3–6/10 | null |
| **Positive control:** PTEN CNV loss → PRAS40 (AKT1S1) T246 | Stouffer z = 3.65 | 10/10 | **passes**. GSK3B S9 z = 2.5 (8/10). The AKT-motif site set is null, so set-level motif tests are diffuse |

**Why the PPP2R2A-protein hit is not a phosphatase effect:**
1. Against the strict null (1,000 random proteins, the same protein in every cohort) it ranks at the 99.5th percentile (empirical p = 0.009). That does not survive the roughly 240 exposure × feature tests.
2. The random proteins that score higher than PPP2R2A are actin (P60709), tubulin (P07437) and other cytoskeletal or housekeeping proteins.
3. GO enrichment for proteins carrying the up-sites is cytoskeleton, intermediate filament and cell adhesion. This looks like a cell-state or compartment axis, not B55 activity.
4. The CNV-independent part of PPP2R2A protein carries the whole signal (β 0.0066), while CNV loss itself gives a non-significant 0.0024.
5. pT and +1 Pro, the core B55 features, do not move.

## 3. Kinase-activity bias (core-v3)
The concern was that CDK or MAPK "activities" would rise in PPP2R2A-deleted tumours because dephosphorylation is reduced. **This is not supported.** With PPP2R2A CNV loss, Stouffer z is CDK1 −0.2, CDK2 +0.8, MAPK1 −1.1 and MAPK3 −1.7, all within null.
- HIPK2 (z = −3.3) and DYRK1A (z = −2.7) do fall with PPP2R2A loss, but against a null of other 8p genes they drop to −0.9 and −0.8. This is an 8p-arm effect, not a PPP2R2A effect (fig4).
- MTOR and AKT1 rise (z = 3.3 and 2.4), which partly survives the 8p null (1.8 and 1.6). Treat this as a hypothesis only.
- PPP2R1A hotspot → GSK3B activity is lower (permutation z = −2.9, UCEC, n = 8). The direction matches Haesen 2016 (inhibitory S9 up), but GSK3B S9 itself is null (z ≈ 0).

## Negatives and caveats
- There is no phosphatase perturbation data that allows study-level held-out learning for PP1/PP2A. Data from Kruse 2020, Cundell 2016 or Kauko 2020 would need to be added.
- Genetic exposures are mostly shallow deletions (median effect about −0.4 log2 protein). PPP2R1A hotspots number 8 in one cohort.
- Site-level meta-p values are anti-conservative because sites are correlated. Rely on the gene and protein nulls instead.
- Requiring a site to be observed in ≥50% of tumours drops most DNA-damage-response sites.
- Unexpected hit: PPM1D gain → pY sites up (β 0.0017, FDR 0.03, 9/10 cohorts). There is no mechanistic expectation for it, so treat it as possible amplicon confounding (17q23).

**Self-score:** Compelling 2/5 · Novel 4/5 (the tumour natural-experiment test of phosphatase rules has not been done) · Promising 2/5. The most useful outcome is the negative control for the field: phosphatase lesions in tumours do not bias kinase-activity scores. The best next step is to add Kruse 2020 B55/B56-depletion data for a proper rule model, then re-test PPP2R2A deep deletions in larger cohorts (CPTAC-3 prostate, METABRIC-matched RPPA).
