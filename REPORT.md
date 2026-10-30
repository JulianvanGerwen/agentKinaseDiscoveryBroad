# REPORT — Discoveries in human kinase signalling from phosphoproteomics (discoveryBroad2)

*AI-run project (Claude Opus 5.5 PI with Opus/Sonnet sub-agents), 2026-10-29. Every claim below quotes the wording of an independent audit, not the producing branch's self-score. No result has been verified by a human. The research record is `RESEARCH_TREE.md`; process rules are in `LESSONS.md`; the literature base is `LITERATURE.md`.*

## Data and tools
- **Perturbation compendium:** 1,937 conditions from 178 studies, including 850 kinase-inhibitor conditions, 172 time series and 86 dose series (decryptM, Hijazi).
- **CPTAC tumours:** 10 cohorts, 1,030 tumours and 523 normal adjacent tissues (NATs), with phospho, protein, RNA, CNV, mutation and clinical data. GDC excision-to-freezing times are available for 7 cohorts.
- **Public data for confirmation:** Menyhart 2023 (CK2 chemical genetics), CPTAC LUAD Confirmatory, Chen 2020 Taiwanese LUAD, a Kuster-lab KRAS dose dataset, GDSC, and DepMap/CCLE.
- **Kinase–substrate predictor:** kinaseSubMain release predictions, 503 kinases.
- **Activity scores:** kinaseActivityMain core-v3 plus an independent PSP mean-z score.
- **Resource limit:** at most a quarter of the node (16 threads), stricter than the brief's "half", following the user's standing instruction of 2026-10-28.

## Ranking scheme
- **Scores:** each candidate was scored 1–5 for Compelling (C), Novel (N) and Promising/relevant (P). Scores below come from an independent auditor who recomputed the result with their own code. Where no numeric audit score exists, the table says whose score it is.
- **Tier 1 (discovery):** audited C ≥ 3.5, replicated in independent data (another lab, cohort or public dataset), and literature-checked.
- **Tier 2 (supported finding):** audited and qualified (C ≈ 3), with replication across cohorts or labs, but either weaker, methodological, or with residual caveats.
- **Work done:** 3 infrastructure tasks, 18 discovery branches, 4 targeted validations, 6 independent audits and 3 reviews (two milestone retrospectives and a final fact-check). Most results were negative and are listed in §4.
- **Score correction:** the final fact-check (`docs/retro/FINAL_factcheck_retro.md`) found that an earlier draft of this table overstated several scores. The table below is corrected.

## 1. Findings

| # | Tier | Finding | C / N / P (source) | Evidence base |
|---|---|---|---|---|
| K1 | **1** | KPNA4 S60 (importin-α3) is a CK2α-dependent, fast-responding site; a published "MAPK-dependent" attribution is not supported | 3.5 / 3 / 3.5 (audit_D); confidence then raised to "high" by audit_F's public rescue data | Compendium (3 cell lines, 2 labs) + independent public chemical-genetic rescue |
| K2 | 2 | Signalling programmes learned from cell-line perturbations are coherent and sign-preserved in human tumour and normal tissue; tumours gain DNA-damage/replication-stress coherence | Coherence: audit_I10 "qualified" (no numeric score). Q8 tissue comparison: unaudited, self-score 4 / 3 / 3 (DNA-damage gain N2) | 10 tumour cohorts, 6 NAT cohorts, measured ischaemia times (MAPK axis) |
| K3 | 2 | Primary-direction responses replicate across labs and cell lines more strongly than reverse-direction responses, which replicate weakly and mainly within a study's culture | 3 / 4 / 3 (audit_I19) | audit_I19 mirror test; consistent with audit_I10's classifier test and the tumour PI3K test |
| K4 | 2 | In tumours, receptor pY that tracks receptor abundance does not reach ERK, whereas acute receptor perturbations in cell lines do; genomic lesions shift cascade offsets, not gains | audit_I24: "qualified, reworded". The 3 / 3 / 4 is the PI's score, not the auditor's | 7 tumour cohorts; cell-line contrast; powered negative |

### K1. KPNA4 S60 is a CK2 site, not a MAPK site
**Evidence chain**
- **Nomination (I03):** integrated perturbation evidence ranked KPNA4 S60 as a CK2 site.
- **Compendium:** the site falls under CX-4945 (z −2.5 to −7.0 in HeLa, HL60 and NTERA2; rank 63 of 12,455 sites in mitotic HeLa, where CK2 inhibition is verified). It falls in 0 of 25 ERK-inhibited experiments.
- **Motif:** its sequence (…ED**S**DID…) matches the canonical acidic CK2 motif. CX-4945 off-target motifs are excluded.
- **Fast-turnover subset:** KPNA4 S60 belongs to a 14-site set that responds within 1 h (64% have a strong CK2 motif; OR 16, p = 1.5e-6), unlike most known CK2 sites, which do not move at 1 h.
- **Independent confirmation, public data (Menyhart 2023, PXD038050; U2OS):**
  - The selective CK2 inhibitor SGC-CK2-1 lowers S60 by −3.73 log2 at 4 h (FDR 3e-4).
  - Expressing an inhibitor-resistant CK2α restores 80% of it (+2.98, FDR 4e-4).
  - All 6 measured members of the fast set behave the same way.

**Relevance.** Hazawa et al. 2019 (*Oncogene*) is the only functional paper on this site. It called the phosphorylation MAPK-dependent, based on a 24 h U0126 experiment and a general phospho-serine pull-down. Our acute data and the chemical-genetic rescue place it under CK2α, so the 24 h MEK effect is probably indirect. For questions about KPNA4 S60 function, the relevant kinase to perturb is CK2, not MEK. We did not check the clinical status of CK2 inhibitors.

**Open questions and caveats**
- Direct phosphorylation by CK2 is likely from the motif but not shown.
- Kinetics within 1 h are shown only with CX-4945.
- Only cell lines were tested.
- Menyhart 2023 found that most sites lowered by CX-4945 are not truly CK2-dependent. This is why the rescue data, not CX-4945, carry the claim.

**Files:** `results/val_I03/`, `results/audit_D/AUDIT.md`, `results/audit_F/AUDIT.md`. Key figure: `results/val_I03/figures/fig_T2_kpna4.png`.

### K2. Perturbation-learned signalling wiring is preserved in human tissue
**Findings**
- **Consensus footprints transfer between studies:** 49 cell-line perturbation footprints, validated on held-out studies (median AUROC 0.71; EGF 0.96, DNA damage 0.97, MEK inhibitors 0.88).
- **Coherent in all 10 CPTAC tumour cohorts.** This survives:
  - protein adjustment;
  - removal of co-abundant and complex site pairs;
  - strict matched nulls (z ≈ 20–30).

  Removing the top 5 global PCs roughly halves the MAPK-axis coherence (MEKi 0.153 → 0.068), but it stays significant in 10/10 cohorts. Coherence is about 4× that of a typical PhosphoSitePlus kinase–substrate set.
- **Already present in normal tissue (Q8, not independently audited):** the MAPK/RTK axis is PC1 in every NAT cohort. 529 signed relationships derived from the compendium agree in sign 91% of the time in both NATs and tumours.
- **Tumour-specific change:** the only robust one is a *gain* in coherence of the DNA-damage/replication-stress programme (meta z 5.3–6.3). It survives composition and proliferation adjustment, and adjustment for the ischaemia *proxy* score. It was not re-tested against measured ischaemia time, and the DNA-damage footprint correlates weakly with measured time (ρ +0.14 [0.07, 0.21]). It newly couples NUMA1 S395 and TP53BP1 S580.

**Pre-analytical artefact (MAPK axis)**
- An audit first suggested the MAPK axis might be pre-analytical.
- Measured excision-to-freezing times (GDC; 7 cohorts, n = 723) are uncorrelated with the MEK/EGF/TNF footprint scores (pooled |ρ| ≤ 0.09). Adjusting for them leaves coherence unchanged.
- The proxy "ischaemia signature" used in that audit tracks the MAPK/stress axis, not measured time.

**Relevance.** Cell-line perturbation experiments can serve as a reference map for interpreting tumour and tissue phosphoproteomes. What distinguishes tumours is a more coherent DNA-damage-response state, not rewired growth-factor signalling.

**Limits**
- The footprints do not beat kinase scores for predicting genotype.
- Three cohorts lack ischaemia times.
- Warm ischaemia is untested.

**Files:** `results/b_I10/`, `results/audit_I10/`, `results/val_B/`, `results/val_E/`, `results/b_Q8/`. Key figures: `results/b_Q8/figures/F1_coherence_NAT_vs_tumour.png` and `results/b_Q8/figures/F2_edge_preservation.png`.

### K3. Primary responses replicate across labs more strongly than reverse responses
**Findings**
- **Primary direction replicates across labs, and more so with effect size.** Primary means decreases under kinase inhibitors and increases under stimuli. Reverse-direction changes also replicate above base rate, but weakly and with a flat dependence on effect size.
- **Mirror test:** across 52 stimulation studies, increases replicate (log odds 0.56 → 1.62 with effect size). Decreases replicate only weakly and flatly (0.48 → 0.59; difference at large effects −0.70 [−1.21, −0.25]). The inhibitor arm shows the reverse pattern, though it leans on MEK inhibitors.
- **Consistent results from other analyses:**
  - As a classifier: down-site footprints transfer better than up-site footprints for inhibitors, and up-site footprints for stimuli, in 75% of footprints each.
  - In tumours: only PI3K-axis down-sites predict PIK3CA/PTEN status (AUROC 0.56–0.57).
- **The reverse layer does not carry across labs, even for the same cell line.** MCF7 measured by two labs disagrees. Within one study (Hijazi), reverse responses do group by cell line, and baseline RNA, protein and DepMap profiles do not predict them.

**Explanations tested and not supported:** timing, cell-cycle redistribution, energy stress, Hippo/density, normalisation, and drug resistance. The cell-cycle, culture-state and resistance tests were pre-registered; the resistance test was underpowered. Calling the reverse layer "compensatory" is an interpretation, not a demonstrated mechanism.

**Relevance.** This is a methods-level caution. Compensatory or feedback increases reported from a single phosphoproteomics study often fail to transfer to another lab, even in the same cell line. Footprints and libraries should be built from primary-direction changes.

**Files:** `results/b_I19/`, `results/audit_I19/`, `results/b_I10/`, `results/audit_I10/`, `results/b_I05/`, `results/b_I22/`, `results/b_I25/`. Key figure: `results/audit_I19/fig_mirror_inhibitor_vs_stimulation.png`.

### K4. Chronic receptor abundance is buffered; lesions shift offsets, not gains
**Findings**
- **Receptor abundance lowers ERK at fixed receptor pY:** more EGFR/ERBB2 protein goes with lower ERK1/2 activation-loop phosphorylation, −0.15 SD [−0.25, −0.05], p = 0.002, 6 of 7 cohorts. It is unchanged after adjusting for the ischaemia proxy score (−0.16); it was not re-tested with measured ischaemia times, though in V-B measured time was uncorrelated with ERK pTpY (ρ −0.05).
- **Once abundance is removed, transmission is similar:** RTK-driven adaptor pY transmits like RTK-independent adaptor pY (0.31 vs 0.44, p = 0.16).
- **Acute perturbations propagate:** in cell lines, EGF stimulation or EGFR/ERBB2 inhibition passes transitively from RTK to ERK (0.58 observed vs 0.51 expected).
- **Classical measurement noise does not explain it:** a simulated fully transitive cascade with the observed reliabilities predicts the opposite gap (P = 0.01).
- **Powered negative on cascade gains:** 0 of 162 pre-declared tests find a genomic lesion that changes cascade gain (slope). Lesions tested: DUSP4/5/6, SPRY, SPRED, NF1, RASA1, PTEN, INPP4B, TSC, STK11, PPP2R2A, and KRAS/BRAF/EGFR/PIK3CA mutations. Power was 84–94% to detect a 35–40% gain change.
- **Lesions shift the offset instead,** as expected: KRAS +0.6 SD of ERK at a given RTK input; PTEN loss and TSC1/2.

**Relevance.** This fits steady-state ERK feedback (Fritsche-Guenther 2011). It may help explain why receptor amplification is a weak marker of downstream pathway output.

**Files:** `results/b_I24/`, `results/audit_I24/`, `results/val_E/`. Key figure: `results/audit_I24/figures/AF1_decomposition_specifications.png`.

## 2. Caveated observations
These are worth reporting, but fall below the key-finding bar.
- **Candidate mTORC1 sites on the lysosomal docking machinery** (auditor C 2.5–3, N 4):
  - LAMTOR1 S97/S98 falls at the 0.11th percentile under catalytic mTOR inhibitors. Its dose potency matches PRAS40 S183, and it is not affected by S6K inhibitors, although S6K inhibition was not verified in those experiments. Both supporting studies come from one lab.
  - RRAGC S15 has a perfect mTOR motif and falls in 6 of 6 organoid lines. It lies next to the known mTORC1 feedback site S21.
  - Downstream AGC kinases are not excluded.
  - Also from this analysis: EPN1 S435 depends on the ERK pathway (direct phosphorylation untested). A literature-only edge, mTOR→AMPKα1 S347, was rediscovered by the method. (`results/b_I03/`, `results/val_I03/`, `results/audit_D/`)
- **Phosphosite turnover atlas** (auditor 3 / 2.5 / 2.5):
  - How fast a substrate responds to verified acute inhibition of its kinase is a reproducible property of the site (split by lab, 0.27 vs 0.15 for random sites). It is specific to the kinase's own inhibitors and holds for ERK and mTOR, but not CDK.
  - Fast sites are better reporters in tumours.
  - Building the library from fast edges gives a small gain in activity inference (+0.035 when whole labs are left out). That gain holds only in the 16–18 "edited" studies, a subset not declared in advance. Across all studies it is +0.010 and not significant, and a few studies carry it. (`results/b_I28/`, `results/audit_H/`)
- **Cellular potency of ATP-competitive drugs fits Cheng–Prusoff:** downstream phosphosite pEC50 vs ATP-free Kd implies about 2 mM intracellular ATP (1.5–3.6 mM). This rests on n = 5 drugs and is a consistency check, not a measurement. Allosteric MEK inhibitors behave differently, and this replicated in an independent KRAS dataset. Dasatinib engages p38α/MK2 in cells at about 1 µM. (`results/b_I17/`; decryptM 2.0 with 133 drugs was access-restricted)
- **Perturbation audit of PhosphoSitePlus** (auditor 2 / 2 / 2.5):
  - Only 7% of edges are testable; 17% are supported vs 5% for random sites.
  - AKT-annotated basophilic sites that depend on MEK/ERK beat chance (p = 0.002), but RSK and AKT cannot be told apart. New pathway candidates are DBNL S269 and AKAP11 T1100.
  - MAP2K2 T394 is contradicted. (`results/b_I26/`, `results/audit_F/`)
- **Kinase activity in tumours tracks activating-site phosphorylation, not kinase abundance** (81 of 117 kinases). This is partly known (Arshad 2019). Mutations and CNV explain about none of the decoupled activity. (`results/b_I08/`)
- **Cohort-specific:** EGFR-mutant LUAD shows low ERK phosphorylation but high ETV4/5 transcriptional output, robust in CPTAC discovery.
  - The low-pERK half did **not** replicate in CPTAC Confirmatory (−0.03, adequately powered, single-pY peptide only).
  - The Taiwanese cohort is underpowered rather than a failed replication: pERK +0.17 with wide CI, while ETV4 +0.55 is in the predicted direction.
  - (`results/b_I29/`)

## 3. Corroboration and relevance to others
- **Literature:** every key finding was checked against `docs/lit/novelty_checklist.md` (97 rows) and by logged web searches. Closest prior work:
  - Hill 2017 (context-specific edges);
  - Arshad 2019;
  - the P100 CMap and PTMsigDB;
  - CoPheeMap and eKPI;
  - Fritsche-Guenther 2011;
  - Hazawa 2019 (contradicted);
  - Menyhart 2023, the source of the confirming data.
- **Honest overall verdict:** one Tier-1 named discovery (K1) and three Tier-2 supported findings (K2–K4). K3 is methodological; K2's tissue comparison is unaudited.
- **What others can use:**
  - K1 corrects a kinase attribution with therapeutic implications.
  - K2 supports using perturbation compendia as references for tissue phosphoproteomics.
  - K3 tells experimentalists which responses will transfer between labs.
  - K4 and the offset-not-gain negative constrain models of oncogenic signalling.

## 4. Negative and abandoned branches (all documented in `RESEARCH_TREE.md`)
- cis-CNV instrumental-variable tests of kinase→substrate edges (IV ≈ naive ≈ null except EGFR);
- dark-kinase function, where 66% of activity is a family echo;
- natural experiments with phosphatase lesions;
- allele-specific signalling predictions (underpowered in CPTAC);
- reverse-layer size vs drug resistance (pre-registered, wrong sign);
- reverse-layer mechanisms;
- MET/ALK→Hippo; PAK4 S148; most tumour reporter sites;
- cross-drug off-target assignment (already covered by decryptM 2.0);
- kinase "reach";
- immune–PTP mimicry;
- ICA modules unexplained by known kinases.

## 5. Next steps
- **K1:** in vitro CK2α kinase assay on a KPNA4 S60 peptide, plus a ≤1 h SGC-CK2-1 time course with the inhibitor-resistant CK2α rescue.
- **K2:** independent audit of Q8, and re-test of the DNA-damage gain with measured ischaemia times.
- **mTORC1 candidates:** mTORC1 assays on Raptor-KO vs S6K1/2-DKO cells for RRAGC S15 and LAMTOR1 S97/S98.
- **decryptM 2.0 (133 drugs):** access would extend the ATP-competition test.

## 6. Limitations
- Kinase activity inference is family-level.
- CDK, CK2, PLK, ATM/ATR and p38 edges are poorly testable with ≤1 h asynchronous data.
- Two labs (Cutillas, Kuster) dominate the inhibitor data, so replication is counted in labs.
- Tumour analyses are cross-sectional.
- Further public data would raise confidence: decryptM 2.0, CK2 knockout phosphoproteomes, and more LUAD cohorts.
