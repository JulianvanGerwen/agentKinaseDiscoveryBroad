# Audit of I10 candidates B (footprint coherence / MEK–EGF axis) and C (immune → pervanadate mimicry)

Auditor: independent re-implementation. I did not modify any b_I10 files. Scripts are in `scripts/audit_I10/00–10`. Tables are `A1–A9`, and the figure is `fig_audit.png`.

## What I re-implemented
- **Data**: own code on `input_tc_padj` (phospho residualised on own protein), tumours only, sites with ≥25% coverage.
- **Correlations**: site-wise z over observed values, then pairwise-complete r (n≥20). I did not use b_I10's rank-normal matrix with NaN→0 imputation; that imputation shrinks every r, so my absolute values run higher than b_I10's.
- **Footprint sets**: the b_I10 T1 sets. PSP sets come from the compendium site table.
- **Matched null B (stricter than b_I10's)**: random sites drawn from the same bin of coverage × tumour SD × communality (R² on the top-5 PCs of all sites) × residue (pY vs pS/T) × proline+1. 200 draws, footprint signs kept. I also kept b_I10's coverage × SD null (A).
- **Further reference sets**:
  - 150 random single compendium experiments (their top-200 |z| sites, as a "biology-agnostic responsive-site" set);
  - the compendium's cold-ischaemia experiments (Mertins 2014, PMID 24719451; 3 tumour-tissue experiments), turned into a 200-site ischaemia footprint.

---
## Candidate B

### B1. Coherence numbers (A1)
Median over 10 cohorts of signed mean cross-protein r. In parentheses: cohorts with p<0.05 vs null B.

| set | base | co-abundant (|r_prot|>0.3) and STRING-physical pairs removed | top-1 / 3 / 5 PCs regressed | PC5 + pairs removed | null B |
|---|---|---|---|---|---|
| MEKi | .153 (10/10) | .151 | .132 / .075 / .068 | .065 (10/10, z≈19) | .006 |
| EGF | .141 (10/10) | .138 | .110 / .063 / .053 | .048 (10/10) | .007 |
| GF_all | .149 | .145 | .111 / .062 / .046 | .044 (10/10) | .008 |
| TNF | .066 (10/10) | .063 | .050 / .039 / .028 | .029 (10/10) | .002 |
| DNA damage | .064 (10/10) | .063 | .069 / .063 / .059 | .061 (10/10) | .006 |
| pervanadate | .057 | .056 | .041 / .027 / .019 | .019 (10/10) | .004 |
| metformin (12 cell lines, single study) | .149 | .149 | .157 / .134 / .112 | .115 | .009 |
| PSP MAPK1 / CDK1 / AKT1 | .030 / .037 / .050 | .030 / .035 / .049 | — / — / .017 / .028 / .025 | .016 / .024 / .024 | .015 / .015 / .006 |
| all 66 PSP sets (median, IQR) | .034 (.018–.055) | | | .017 | |

- **The numbers reproduce.** The sets are coherent in 10/10 cohorts. They stay so under every confound I tried, including the null matched on communality, residue and proline+1, which b_I10's null did not match.
- **Same-protein and same-complex pairs are not the cause.** Same-protein pairs were already excluded. Removing pairs whose proteins co-vary (|r|>0.3) or bind physically (STRING ≥700) changes coherence by ≤3%.
- **Selection bias (2c) is not the cause.** Footprint sites match random sites on coverage (~0.72), SD and communality (~0.29). Re-matching on variance, missingness, residue, P+1 and communality leaves z ≈ 20–30.
- **Global factors (2b) are the big effect.** Regressing out the top 3–5 all-site PCs removes about 55–65% of the MAPK/EGF/TNF/pervanadate coherence: MEKi drops from .153 to .068 and EGF from .141 to .053. DNA damage and metformin are essentially unaffected.
  - So the MAPK-type footprints sit largely on the dominant tumour phospho PCs. Their own scores correlate |ρ| 0.33–0.60 with all-site PC1/PC2 (A9).
- **Specificity (A1 random experiments).** Top-200 sites of a random compendium experiment have median coherence .008. The MEK, EGF, GF and metformin footprints beat 100% of them, and TNF and DNA damage beat ≈95%.
  - **The best random experiment was cold ischaemia (.127).** That is what prompted the ischaemia test in B3.
- **PSP comparison (2d).** It is fair on coverage, SD and communality, and mean pairwise r does not depend on set size. But "≈0.02" picks weak sets: the PSP median is .034, and RPS6KA1, MAPKAPK2 and SRC are .09–.10. A defensible ratio is **≈4× the median PSP set**, both before (.153 vs .034) and after PC5 (.068 vs .017). It is not 5× vs "PSP ≈ .02".
- **Direction-shuffle null.** It is a weak test for MEK and EGF, which are 93% and 90% one sign. Only the mixed-sign sets (TNF, DNA damage) gain real support from it.

### B2. The "−0.88" MEK/EGF correlation (A2)
- **Reproduced.** Own scores give median r = −0.92 (range −0.85 to −0.93).
- **The sign is built in.** The MEK footprint is 186/200 down sites (ERK targets that lose phosphorylation) and EGF is 179/200 up sites (largely the same kind of targets gaining it). With MEK-score = Σ sign·z, the MEK score is approximately the negative of MAPK output, so −0.88 means "two MAPK-target sets correlate +0.88".
  - Even random matched sites with the MEK vs EGF sign patterns give r ≈ −0.34 (median; 5–95%: −0.6 to −0.1).
- **The magnitude is partly shared sites, partly real.** Dropping the 35 shared sites gives −0.85, dropping shared proteins gives −0.84, and MEK-down-only vs EGF-up-only gives −0.85. For comparison, two disjoint halves of the PSP MAPK1/3 substrates correlate only +0.68.
- **Biologically, it is the MAPK/ERK output axis.**
  - MAPK output (−MEK score) vs core-v3 MAPK1/3: ρ 0.83 (tc_padj; range .71–.91).
  - vs PROGENy MAPK (RNA): ρ 0.18.
  - RAS/RAF hotspot AUROC: 0.62 per-cohort median and 0.57 pooled, vs 0.67 / 0.57 for MAPK1/3 kinase activity. It adds nothing over the kinase score, as b_I10 itself found.
- **"PC1 27%" is real but means less than it sounds (A6).**
  - b_I10's 49-footprint PCA reproduces at 27.1%. Matched random sets with the same weights and signs give 3.7% (3.4–4.1%), so PC1 is not a scoring artefact.
  - About 8 footprints are redundant MAPK readouts (EGF, GF_all, EGFR, MET/ALK, MAPKaxis, RAF, p38, MEK). Collapsing them to one still leaves PC1 at 22%.
  - PC1 loadings correlate ρ = −0.77 with the fraction of down sites in each footprint. Up-dominated stimulus footprints (EGF, GPCR, TNF, pervanadate, immune receptor) load +; inhibitor footprints load −.
  - So PC1 is a generic "stimulus-responsive phosphorylation high vs low" axis whose best single readout is MAPK1/3 activity (ρ 0.70).

### B3. Confounders of the axis (A5, A8, A9)
- **Purity, stroma, immune and proliferation are not the explanation.** Medians for MEK/EGF scores are |ρ| ≤ 0.17 with purity, ESTIMATE stromal/immune and E2F/G2M.
- **Cold ischaemia is a major confounder, and "Not ischaemia-driven" in SUMMARY.md is contradicted.** I found no script in b_I10 that tests ischaemia, and CPTAC has no ischaemia times in our assets.
  - The compendium cold-ischaemia footprint shares only 4–8 sites with MEK/EGF. Even so, its tumour score correlates with MEK −0.70, EGF +0.69, TNF +0.76, pervanadate +0.72 and GPCR +0.50 (median of 10 cohorts, all cohorts the same sign). DNA damage (−0.11) and metformin (0.03) do not track it.
  - A second ischaemia proxy (immediate-early/stress phosphosites, unadjusted) gives ρ 0.31 with MAPK output.
  - Regressing the ischaemia score out of every site cuts coherence: MEK .153→.084, EGF .141→.074, TNF .066→.025, pervanadate .057→.024. DNA damage is unchanged.
- **The genotype link survives ischaemia adjustment.** After residualising on the ischaemia score, MAPK output still separates RAS/RAF-hotspot tumours (pooled AUROC .557 vs .568 before). The ischaemia score alone gives .530. Residual MAPK output still correlates with residual MAPK1/3 activity (ρ ≈ 0.72).
- **Interpretation.** About half the tumour MAPK/stimulus-footprint co-variation is shared with a pre-analytical cold-ischaemia response. Ischaemia itself activates ERK/p38/immediate-early signalling, so with these data the two cannot be told apart. What remains after removing it is still coherent and carries the RAS/RAF signal.

### B4. Direction asymmetry (A3)
- **The b_I10 T2 table re-tabulates, but the overall claim depends on footprint mix.** Across all footprints: down .68 vs up .59, with down > up in 70%. In inhibitors: .67 vs .54 (78%). In stimulus/stress footprints: .75 vs .78 (50%).
- **Own LOSO** (independent code, ≥2 training studies, 28 footprints):
  - inhibitors: down .65 vs up .55, down > up in 75%; unchanged when size-matched;
  - stimuli: up .78 vs down .55, **up > down in 75%**.
  - So this is not a down-vs-up asymmetry. **The primary-response direction transfers better** (direct dephosphorylation for inhibitors, induced phosphorylation for stimuli), and the compensatory/secondary direction transfers poorly. That is expected biology and does not replicate "Finding A" as a general down-site property.
- **PI3K in tumours reproduces.** Pooled AUROC: down .564 (p = .01 vs matched random sets), up .460 (two-sided p ≈ .09, i.e. trending inverse), size-matched down .554.
  - Per cohort, the down sites carry it in BRCA, UCEC, LSCC and GBM (.59–.64). COAD is .47 and HNSCC has only 8 positives.
  - The effect is real but small. It is again consistent with "primary direction" (the down sites are PI3K/AKT/mTOR targets).

### Verdict B: **QUALIFIED**
- **Survives**:
  - compendium footprints are coherent co-varying units in tumours, robust to protein adjustment, removal of co-abundant/complex pairs, a re-matched null (variance, missingness, residue, communality) and top-5-PC removal;
  - they are ~4× more coherent than the median PSP kinase–substrate set and beat ~95–100% of random-experiment site sets.
- **Must be changed**:
  1. "−0.88" is largely MAPK1/3 activity measured twice with opposite sign conventions. Report it as "two independent MAPK-target sets correlate +0.85 after removing shared sites", and say the axis ≈ MAPK1/3 kinase activity (ρ 0.83). It does not add to kinase scores for RAS/RAF.
  2. "Dominant tumour axis, not ischaemia-driven" is not supported. The axis correlates |r| ≈ 0.7 with a compendium cold-ischaemia footprint, and ~half its coherence is removable by that score. It should be presented as "MAPK/stress-responsive axis, partly confounded by pre-analytical ischaemia".
  3. "5× vs ≈0.02" should become "≈4× the median PSP set". Some PSP sets (MAPKAPK2, RPS6KA1, SRC) reach the footprint range.
  4. The direction asymmetry is "primary-response direction > secondary", not "down > up".

**Defensible wording**: "Consensus perturbation footprints derived from cell-line phosphoproteomics remain coherent in all 10 CPTAC cohorts. This holds after protein adjustment, removal of co-abundant/complex protein pairs, stringent matched-site nulls and removal of the top 5 global PCs. Signed pairwise r is ~4× that of a typical PSP kinase–substrate set. The strongest units (MEKi/EGF/GF) form a single MAPK-output axis that tracks MAPK1/3 kinase activity (ρ≈0.8) and RAS/RAF status. Roughly half of its tumour-to-tumour variance is shared with a cold-ischaemia response signature and may be pre-analytical."

---
## Candidate C: immune infiltration → pervanadate mimicry (A4, A7)

### Reproduction
With own scoring and OLS on immune + purity, the Stouffer meta-z is 5.4. Five cohorts are positive at p < .05; CCRCC is opposite. b_I10 reported 4.3–4.5 with its nz scores.
- **But generic sites carry a large part of it.** Matched random sets with the same signs give a median meta-z of 3.0 (95th percentile 4.5). The footprint's excess over that null is p = .01.

### The leukocyte pY confound
It is not the mechanism, but only because CPTAC IMAC barely measures pY.
- Only 92–331 pY sites per cohort pass coverage, and only 5–14 of the footprint's 104 pY sites are measured.
- **pY-only set**: meta-z 2.9, computable in just 4 cohorts.
- **pS/T-only set**: meta-z 5.0.
- **Removing HPA immune-cell-specific proteins** (1,708 genes): z 5.4, unchanged.
- **Removing proteins whose abundance correlates with ESTIMATE immune** (4,762, a very liberal cut): z 4.7.
- **Removing both**: z 4.1 (111 sites, 3 cohorts significant, p = .025 vs matched null).

**Drivers (per-site meta-z)** are mainly pS/T sites on non-immune proteins: MARK2 S409, AKAP13 S2563, RAB7A S72, AHNAK S1068/S4850, PNN S690, ITSN2 Y968, HSPB1 S15, RBM34 S14, NOLC1 S563, and STAT3 Y705 (the IL-6/inflammation marker). Classic leukocyte pY (LCK, ZAP70, SYK, LCP2) is absent or unmeasured.

### The real problem is specificity
- In tumours the "pervanadate" score correlates r = 0.85 with the EGF footprint score, −0.80 with the MEK footprint score, 0.70 with TNF and 0.72 with the cold-ischaemia score. Despite only 7–13 shared sites, it is another readout of the PC1 "stimulus/stress-responsive phosphorylation" axis.
- The immune association survives adjustment for MEK + TNF + EGF scores (z 4.3), for the ischaemia score (z 4.7) and for PC1–3 (z 6.1). So immune-high tumours do have higher phosphorylation at these sites beyond the MAPK axis.
- Nothing ties it to PTP inhibition: the pY component is unmeasured and pS/T sites carry the signal. The TNF (z ≈ 3.8) and immune-receptor footprints show the same direction.
- **"Not seen with LCK/ZAP70/SYK/JAK1 kinase scores" is a weak contrast.** Those are protein-adjusted (tc_padj) z_core scores. Infiltration raises immune kinase protein, which the adjustment removes, and their substrates are mostly pY sites that CPTAC barely measures.

### Verdict C: **FAILS as stated (PTP-inhibitor mimicry); a weaker, re-labelled association survives**
**Defensible wording**: "In 5/10 cohorts, immune-infiltrated tumours show higher phosphorylation at a set of pS/T sites induced by phosphatase inhibition in cell lines (HSPB1 S15, AHNAK, AKAP13, MARK2, STAT3 Y705). This holds after removing immune-cell-specific and immune-co-abundant proteins (meta-z 4.1) and after adjusting for MAPK, TNF and ischaemia scores. The tyrosine component cannot be assessed in CPTAC IMAC data, so no PTP-inhibition mechanism is implied." Self-rating C2/N2/P2.

---
## Single most convincing figure
`results/audit_I10/fig_audit.png`
- **Panel A**: median coherence for footprints vs PSP sets in three versions (base; co-abundant/STRING pairs removed + top-5 PCs regressed; cold-ischaemia score regressed), with the strict matched null. Footprint coherence survives everything but halves for the MAPK, TNF and pervanadate units, and DNA damage is untouched.
- **Panel B**: per-cohort correlation of each footprint score with the compendium cold-ischaemia footprint score. MEK/EGF/TNF/pervanadate/GPCR sit at |r| 0.5–0.75 in every cohort; DNA damage and metformin are near 0.

## Recommended follow-ups
1. Get CPTAC ischaemia times (CDAP/PDC biospecimen metadata) and test the axis directly.
2. Report DNA-damage coherence as the cleanest (ischaemia- and PC-independent) example.
3. Treat the metformin footprint (single study, 12 cell lines) with caution: its coherence is high but it is not one of the claimed footprints.

## Files
- **Scripts** (`scripts/audit_I10/`): 00_prep (HPA immune list, STRING), common, 01_coherence, 02_axis, 03_direction, 04_pervanadate, 05_ischaemia, 06_pc1_null, 07_pervanadate_specificity, 08_mapk_isch, 09_correlates, 10_figure.
- **Results** (`results/audit_I10/`):
  - A1_coherence.tsv, A1_random_experiments.tsv
  - A2_axis_percohort.tsv, A2_correlates.tsv, A2_pca.tsv
  - A3_T2_retab.tsv, A3_loso_own.tsv, A3_pi3k_tumour.tsv
  - A4_pervanadate.tsv, A4_site_drivers.tsv
  - A5_ischaemia.tsv, A6_pc1_null.tsv, A7_pervanadate_specificity.tsv, A8_mapk_vs_ischaemia.tsv, A9_correlates.tsv
  - fig_audit.png
