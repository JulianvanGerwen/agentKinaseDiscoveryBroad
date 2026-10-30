# Audit F: I26 PSP-edge corrections (candidate F) and KPNA4 S60 under selective CK2 inhibition

Auditor: independent re-analysis. I did not modify any b_I26 file. My own code is in `/workspace/scripts/audit_F/` (01–06; the `t2_*` scripts are the Task 2 retrieval). Tables are in `/workspace/results/audit_F/tables/`.

## Bottom line

- **Task 1 (F).** Only part of the claim holds.
  - The **AKT → MEK/ERK-pathway** pattern is real and beyond chance: 10/30 AKT-unsupported units versus 8% for matched random basophilic sites (p≈0.002).
  - The **CDK → ERK part is at chance level**: 5.3% versus a 4.9% null. So is mTOR.
  - Of the 11 named sites:
    - 4 are supported corrections or reassignments: EPHA2 S897, CHEK1 S280, DBNL S269, PIK3C2A S259 (the last only as a context addition).
    - 1 is moderately supported: AKAP11 T1100.
    - 4 are ambiguous: CTTN S405, DNMT1 S143, EDC3 S161, FOXK2 S373.
    - 2 are not supported: **MAP2K2 T394** and **TFPT S180**.
  - Two of the four "literature-consistent" corrections are not new: **PSP already lists RSK for EPHA2 S897 and CHEK1 S280**, and **ERK for CTTN S405**. The fourth, MAP2K2 T394, **contradicts** its primary literature. Banks 2015 shows the site is lost on acute analogue-sensitive Cdk5 inhibition in HEK293 cells. The "ERK feedback" attribution comes only from a vendor datasheet.
  - Neither motifs nor the AKT-inhibitor data can tell RSK from AKT, for two reasons:
    - Cantley percentiles for AKT and RSK are essentially identical on these sites.
    - In 179 of 201 AKT-inhibitor site-conditions where inhibition was verified, MAPK markers *rose* (relief of feedback). A dual AKT/RSK site would therefore look flat.
- **Task 2 (KPNA4 S60).** **Public data found, and the result is decisive.** In Menyhart et al. 2023 (Curr Res Chem Biol 3:100041; U2OS CSNK2A1-WT vs inhibitor-resistant triple mutant):
  - The selective probe SGC-CK2-1 (5 µM) lowers KPNA4 S60 by **−3.73 log2 at 4 h** (FDR 3e-4) and **−3.58 at 24 h**.
  - The inhibitor-resistant CK2α **restores it by +2.98 log2** (H/M, FDR 4e-4). That is about 80% rescue.
  - I checked these numbers directly against the xlsx.
  - So S60 is a CK2α-dependent site. This upgrades audit D. Turnover speed (≤1 h) is still not tested, because the earliest time point is 4 h.

---

## Task 1: candidate F

### Method (independent of I26's PSP-gold control)

- **Conditions.** Single-agent kinase-inhibitor conditions with QC ok, not combinations, and ≤8 primary targets. They are classed by primary target: AKT, PKC, CDK1/2/5/9, mTOR, MEK, ERK, RAF, RSK, ALK/ROS1, EGFR. Values are robust z per condition, with the sign corrected.
- **Verified inhibition.** I used hand-picked canonical readouts, which are not the PSP-gold sets I26 used. The Stouffer score had to be ≤ −2.5:
  - AKT: GSK3A/B S21/S9, PRAS40 T246, TSC2 T1462, FOXO3 S253, FOXO1 S256, AS160 T642.
  - MAPK: ERK1/2 activation loop, plus RAF1/BRAF/SOS1/EGFR feedback sites.
  - CDK: RB1 S780/795/807/811, NPM1 T199, LMNA S22.
  - PKC: MARCKS, PKD1 S738/742.
  - mTOR: 4E-BP1, ULK1 S758, PRAS40 S183, RPS6 S240/244.
- **Labs.** Studies were mapped to last-author labs using PubMed esummary.
  - Hijazi 2020 and Wilkes 2015 are both from the **Cutillas** lab.
  - Bekker-Jensen 2020 and Emdal 2018 are both from **Olsen**.
  - decryptM, Koch 2016, Ruprecht 2017 and Reinecke 2019 are all from **Kuster**.
  - Basken 2018 and Stuart 2015 are both from **Ahn**.
  - Consequence: I26's "≥2 studies" is often only 2 labs. Verified AKT inhibition comes from **only 3 labs** (Cutillas, Kuster, Olsen). Verified CDK inhibition comes from **1 lab** (Burgess: RO3306 at mitotic exit). RSK inhibitor conditions (BI-D1870, Roux) could not be verified by markers.
- **Excluded study.** Thibault 2017 (vemurafenib in KRAS-mutant HCT116, a paradoxical-RAF-activation design) is removed from MAPK evidence. In that study the test sites mostly *rise*.

### (a) Per-site lab counts under verified inhibition

`tables/site_lab_counts_final.tsv`. Each cell gives labs down (lab mean z ≤ −1) / labs flat (> −0.5) / labs tested.

| Site | PSP group | Annotated-kinase inhibitors | MEK/RAF | ERK inhibitors | Other | Leave-one-lab-out |
|---|---|---|---|---|---|---|
| EPHA2 S897 | AKT | AKT **0/3/3** (+0.8, +1.1, +2.4) | **4/1/5** | 2/0/2 | EGFR 3/0/3 | Robust: ≥3 labs down whichever lab is dropped; AKT flat in all 3 |
| CHEK1 S280 | AKT | AKT 0/1/2 (Cutillas −0.6, Kuster +1.0) | 3/2/5 | **3/0/3** | Roux: RSK1/2 RNAi −3.75 and BI-D1870 −3.25 (HEK293+PMA) | MAPK side robust. The "not AKT" side rests on Kuster alone; without it, Cutillas −0.6 is indeterminate |
| DBNL S269 | AKT | AKT **0/3/3** (−0.2, +2.3, +3.1) | **4/1/6** | 2/1/3 | ALK 2/0/2 | Robust (≥3 labs down whichever lab is dropped) |
| PIK3C2A S259 | CDK1 | CDK 0/1/1 (Burgess only) | **6/0/7** | 2/0/3 | ALK 2/0/3 | MAPK side robust. CDK1 is effectively untested (1 lab, mitotic-exit design) |
| AKAP11 T1100 | mTOR | mTOR **0/4/4** | 3/2/6 | **3/0/3** | – | ERK-inhibitor side survives (2/2 after any drop). MEK/RAF side is weaker |
| CTTN S405 | PKCδ, CDK1/2 | PKC 0/2/2; CDK 0/1/1 | 3/1/4 (two labs at −8) | **1/2/3** (Cutillas +4.0, Kuster +1.5) | – | ERK-inhibitor data contradict ERK. Only MEK inhibitors support it |
| DNMT1 S143 | AKT | AKT **1/1/2** (Cutillas −1.0, Kuster +0.4) | 3/2/5 | 2/0/3 | Roux RNAi −1.25, BI-D1870 ≈0 | Not robust: AKT evidence is split |
| EDC3 S161 | AKT | AKT **1/1/3** (Cutillas −1.1, Olsen −0.8, Kuster +1.6) | 3/2/6 | 1/1/2 | ALK 3/1/4 | Not robust. ALK inhibitors hit both PI3K/AKT and MAPK |
| FOXK2 S373 | CDK1/2 | CDK 0/1/1 (Burgess +1.8) | 0/1/2 | 0/1/1 | mTOR 2/0/2; AKT −2.0 (Cutillas) | mTOR rests on 2 labs. It is an mTOR/AKT-pathway site, not ERK |
| MAP2K2 T394 | CDK5 | CDK 0/1/1 (CDK1 inhibitor only, mitotic exit) | 3/3/7: MEK inhibitors flat in Kuster, Olsen, Cutillas | 1/1/3 | – | Not supported |
| TFPT S180 | CDK1 | CDK 0/1/1 (+3.3) | **0/5/6** | 0/2/3 | – | Not supported. I26's RAF call came from Thibault (−1.1) and Trajanoski (−0.3) |

Notes:
- The I26 per-site calls almost all rest on n=2 studies on the PSP side (`psp_n_meas`=2). For AKT and CDK these are Hijazi (Cutillas) plus decryptM (Kuster), or Burgess.
- For CDK edges there is no valid test of the annotated kinase in interphase cells. The only verified CDK experiment is CDK1 inhibition at mitotic exit. "Unsupported CDK1/5 edge" is therefore not evidence of a wrong edge.

### (b) Can the data distinguish RSK from AKT?

**Motif.** `tables/motif_percentiles.tsv`. On the basophilic sites the AKT and RSK percentiles are indistinguishable:

| Site | AKT (max) | RSK (max) |
|---|---|---|
| CHEK1 | 95 | 95 |
| EPHA2 | 94 | 99 |
| DNMT1 | 94 | 98 |
| EDC3 | 100 | 100 |
| DBNL | 96 | 94 |

Across all sites:
- AKT1 vs RPS6KA1 percentiles have Spearman ρ = 0.82.
- Of sites at ≥90 for either kinase, 54% are ≥90 for both.
- The motif adds no discrimination. The "RSK" in F's wording is a pathway inference: the site is MEK-dependent and motif-compatible with RSK.

**Was AKT inhibition verified where the site did not fall?** Yes. `tables/akt_vs_mapk_verification.tsv` and `akt_verified_conditions_per_site.tsv` show that the median AKT-marker Stouffer score is about −6.8 in AKT-inhibitor conditions. The site did not fall (z > −1) in:

| Site | Verified AKT-inhibitor conditions where the site did not fall |
|---|---|
| CHEK1 | 12/14 |
| EPHA2 | 16/17 |
| DBNL | 19/22 |
| DNMT1 | 9/14 |
| EDC3 | 15/18 |

In MEK/ERK/RAF conditions with verified inhibition, the AKT markers did *not* fall (median ≈ 0). So the MEK-inhibitor effects are not routed through AKT.

**Major confound.** In **179/201** AKT-verified site-conditions, the MAPK markers **rose** (median +2.4): AKT inhibition relieves feedback on ERK/RSK. A site phosphorylated by both AKT and RSK would therefore stay flat under AKT inhibitors. "Flat under AKT inhibitors" cannot exclude AKT.
- In the few AKT-verified conditions with no MAPK rise (n=1–3 per site), DNMT1 (−1.0), EDC3 (−0.9) and DBNL (−0.75) drift down.
- One same-lab, same-cell test exists: decryptM A549 at 30 min, with MK2206, selumetinib, and the combination (`tables/decryptM_A549_combo.tsv`).
  - DBNL S269 rises under MK2206 alone (+1.0 log2); the combination abolishes the rise (0.03). EPHA2 behaves the same way (+0.49 → −0.16).
  - This fits MEK/RSK dependence with AKT-inhibitor-induced feedback.
  - A549 (KRAS-mutant) barely responds to selumetinib at these sites, so the test is weak.
- The only direct RSK perturbation is Roux 2014 (one lab). It gives RNAi and BI-D1870 support for **CHEK1 S280** only. DNMT1 is weak (RNAi −1.25, BI-D1870 ≈ 0). EDC3 and the others are flat or unmeasured.

**Conclusion on (b).** The data can show "MEK-pathway-dependent and AKT-inhibitor-insensitive in acute settings". They cannot show "RSK rather than AKT is the direct kinase", except for CHEK1 S280 (one lab with RSK perturbation, plus the literature). The ALK-inhibitor-based calls (EDC3, and AGO2/HSPB1/CARHSP1 in the full list) cannot separate AKT from RSK at all, because ALK inhibitors suppress both arms.

### (c) Is the pattern beyond chance?

`tables/null_tests.tsv`, script `03_null.py`. The same I26 rule was applied: alternative group "supported", ≤2 groups supporting the site, alternative mean ≤ −1.5, MAPK-axis group, motif ≥85.

| Null | PSP group | Observed | Null mean [5–95%] | p |
|---|---|---|---|---|
| A: matched random sites (AKT-unsupported, not PSP-AKT; matched on residue, AKT-motif bin, n studies) | AKT (30 units) | **0.33** | 0.076 [0, 0.17] | **0.002** |
| A + volatility matching | AKT | 0.33 | 0.085 [0, 0.17] | 0.002 |
| A | CDK1/2/5/9 (132) | 0.053 | 0.049 [0.015, 0.077] | 0.46 |
| A + volatility | CDK | 0.053 | 0.046 | 0.36 |
| A | mTOR (33) | 0.061 | 0.069 | 0.68 |
| B: permuted inhibitor labels (each group's verified conditions replaced by an equal number drawn at random) | all 270 unsupported units | 0.085 (MAPK + motif) | 0.031 [0.015, 0.056] | 0.003 |

How to read this:
- The abundance of MEK/RAF experiments does **not** explain the AKT pattern.
  - PSP-AKT basophilic sites that are flat under AKT inhibitors are about 4× more often MEK-dependent than matched random basophilic sites.
  - This fits PSP AKT edges that are partly RSK sites, from in-vitro AGC promiscuity or historic attributions.
- The **CDK → ERK part of the headline is at chance level**. Given how many MAPK experiments there are, about 5% of any CDK-unsupported SP site gets an "ERK" alternative.
- The permutation null confirms that the overall MAPK-axis excess is non-random (p ≈ 0.003). The AKT units drive it.

### (d) Literature and PSP status

Read from abstracts and full text where reachable. PSP status comes from the local `Kinase_Substrate_Dataset`.

| Site | Claimed | What the paper says | Alternative already in PSP? | Verdict on "literature-consistent" |
|---|---|---|---|---|
| CHEK1 S280 | RSK not AKT | Li et al. 2012 MBoC 23:1582 (PMID 22357623): serum-induced S280 phosphorylation is p90RSK-dependent; RSK phosphorylates S280 stoichiometrically in vitro; "p90 RSK, but not Akt/PKB". Other groups (melanoma, Oncogene; glucose deprivation, Cancer Sci) also show RSK→S280 | **Yes**: RPS6KA1 (in vivo and in vitro). I26's own class column wrongly says "novel alternative" | Consistent, but a **confirmation of an existing PSP edge**, not a discovery |
| EPHA2 S897 | RSK | Zhou et al. 2015 **Nat Commun** 6:7679 (PMID 26158630; I26 wrote "Mol Cell"): cytokine-induced S897 phosphorylation is RSK- and not AKT-dependent. AKT edge from Miao 2009 | **Yes**: RPS6KA1 and RPS6KA3 | Consistent. Already in PSP. AKT is context-specific, not wrong |
| CTTN S405 | ERK | Established ERK site (Campbell 1999; Martinez-Quiles 2004; Kelley 2010 PMC2973953). PRKCD edge is in vivo only; CDK1/2 in vitro only | **Yes**: MAPK1 and MAPK3 | Consistent in the literature. Our data are **mixed**: ERK inhibitors *raise* S405 in Cutillas (+4.0) and Kuster (+1.5) |
| MAP2K2 T394 | ERK feedback, not CDK5 | Banks et al. 2014/2015 Nature 517:391 (PMC4297557): Cdk5/p35 spiked into lysates phosphorylates mouse MEK2 T395/T397 (= human T394/T396); **acute inhibition of analogue-sensitive Cdk5 in HEK293 abolishes pT394**; T→A mutants hyperactivate MEK. "ERK2 feedback at T394" appears only on a vendor datasheet, which looks like a conflation with MEK1 T292 | No. PSP: CDK5 → mouse Map2k2 T395, in vitro | **Contradicted**. The original evidence is a chemical-genetic cell experiment, and the acute data here never tested CDK5 |

The PSP origin of AKAP11 T1100 is mouse in-vivo Mtor → Akap11 T1092. The human acute mTOR-inhibitor data (4 labs flat) argue it is not an acute mTOR-dependent site in human cell lines.

I did not independently verify the "weak original evidence" notes for DNMT1, PIK3C2A, FOXK2 and EDC3 beyond the PSP flags. PSP lists DNMT1 S143 AKT1 as in vivo and in vitro, and PIK3C2A S259 CDK1 as in vivo and in vitro.

### Verdict per site

| Site | Verdict | Defensible statement |
|---|---|---|
| EPHA2 S897 | **Supported** (not novel) | MEK/EGFR-pathway-dependent and AKT-inhibitor-insensitive in 3 labs. This agrees with the existing PSP RSK edge. The AKT edge is context-specific |
| CHEK1 S280 | **Supported** (not novel) | MEK/ERK- and RSK-dependent (Roux RNAi/BI-D1870). This agrees with Li 2012 and the existing PSP RSK1 edge. The "not AKT" side rests on 1–2 labs |
| DBNL S269 | **Supported as pathway reassignment** | MEK-dependent in 4/6 labs; flat or rising under verified AKT inhibition in 3 labs. RSK is a candidate by motif only (AKT motif equally high). The AKT edge is not supported in acute settings |
| PIK3C2A S259 | **Supported as an addition, not a correction** | MEK-dependent in 6/7 labs. The CDK1 (G2/M) edge was not tested in interphase |
| AKAP11 T1100 | **Moderately supported** | ERK inhibitors lower it in 3/3 labs; MEK/RAF 3/6; mTOR inhibitors flat in 4/4. A candidate ERK site (TP, ERK motif 85–88) |
| CTTN S405 | **Ambiguous** | Literature ERK, already in PSP. PKCδ not supported (2 labs). Our ERK-inhibitor data conflict |
| DNMT1 S143 | **Ambiguous** | The AKT inhibitor lowers it in Cutillas. MEK dependence is in 3/5 labs. RSK perturbation is weak |
| EDC3 S161 | **Ambiguous** | The AKT inhibitor lowers it in 2/3 labs. ALK-inhibitor evidence cannot separate AKT from RSK |
| FOXK2 S373 | **Ambiguous** | mTOR/AKT-pathway-dependent (2 labs). The direct kinase is unknown. Not an ERK/RSK case |
| MAP2K2 T394 | **Not supported** | Primary literature (Banks 2015) supports CDK5 with cell-based chemical genetics. MEK inhibitors flat in 3 labs |
| TFPT S180 | **Not supported** | No MAPK-inhibitor dependence (0/6 labs) once the paradoxical-activation study is excluded |

### Overall verdict on F

The headline is overstated:
- "27 likely corrections" counts 2-study calls in which the annotated kinase was often not testable (all CDK edges).
- The CDK half of "AKT/CDK → ERK/RSK" is at chance level.
- Two of the four flagship literature cases are PSP edges that already exist, and one (MAP2K2) is contradicted.

What survives:
- PSP-AKT basophilic sites that do not respond to AKT inhibitors are enriched about 4× over matched random sites for MEK-pathway dependence.
- A handful of sites (EPHA2, CHEK1, DBNL; PIK3C2A and AKAP11 as pathway additions) are robust to leave-one-lab-out.

**Defensible wording.** "In acute inhibitor phosphoproteomics, a subset of PSP AKT-annotated basophilic sites (about 1/3 of testable, AKT-unresponsive units; about 4× matched background) depend on MEK/ERK signalling and not on AKT activity. This is consistent with RSK being the relevant kinase in these contexts (EPHA2 S897 and CHEK1 S280 recapitulate published RSK assignments already in PSP). DBNL S269 and AKAP11 T1100 are new MEK/ERK-pathway-dependent candidates. Motif and acute data cannot distinguish RSK from AKT directly, and CDK edges could not be tested."

**Auditor scores: C 2 / N 2 / P 2.5.**

---

## Task 2: KPNA4 S60 under selective CK2 inhibition

Retrieval was delegated to a sub-agent (about 9 minutes, inside the 40-minute time box). I independently re-read the key KPNA4 rows from the downloaded xlsx with `python -I` (`05_spotcheck_menyhart.py`). The DOI resolves to the same Elsevier PII as the files.

**Dataset.** Menyhart et al. 2023, *Curr Res Chem Biol* 3:100041, doi 10.1016/j.crchbi.2023.100041 (PRIDE PXD038050).
- Cells: U2OS Flp-In T-REx expressing CSNK2A1-WT or the inhibitor-resistant V66A/H160D/I174A mutant (TM).
- Triple SILAC: L = WT + DMSO, M = WT + inhibitor, H = TM + inhibitor.
- Treatments: CX-4945 30 µM or SGC-CK2-1 5 µM, at 4 h and 24 h, n=6.
- Files: `/workspace/cache/audit_F/dl/menyhart2023/mmc4.xlsx` (CX-4945) and `mmc5.xlsx` (SGC-CK2-1). Table: `/workspace/results/audit_F/task2_menyhart2023_module_sites.tsv`.

**KPNA4 S60.** Multiplicity 1, window …HEDICED**S**DIDGDYR…, so the numbering matches directly.

| Arm | M/L (drug effect) | FDR | H/L (resistant + drug) | H/M (rescue) | FDR |
|---|---|---|---|---|---|
| SGC-CK2-1 4 h | **−3.73** | 3e-4 | −0.82 | **+2.98** | **4e-4** |
| SGC-CK2-1 24 h | −3.58 | 0.038 | −0.79 | +2.83 | 0.074 |
| CX-4945 4 h | −3.12 | 2e-4 | −2.12 | +1.05 | 0.028 |
| CX-4945 24 h | −3.18 | 5e-4 | −2.44 | +0.78 | 0.19 |

- **Sanity check.** PSP CK2 substrates fall in every arm (median M/L −0.57 to −1.33). They are rescued only in the SGC arm (H/M +0.78 / +1.16). The CX-4945 arm rescues poorly across the board, which is consistent with CX-4945 off-targets or incomplete resistance of TM to 30 µM. So the experiment worked.
- KPNA4 S60 is missing from the authors' curated CSNK2A1 substrate list only because the residual H/L (−0.82) misses their ±0.585 "Rescued" window.

**The 14-site CX-4945-specific module (audit D).**
- 6 of the 14 sites are quantified, and every one falls under SGC-CK2-1 with a positive rescue (H/M +1.9 to +4.0):
  - KPNA4 S60
  - ANLN S54 (SGC 4 h −2.91, H/M +3.55)
  - FXR1 S409 (−2.7, H/M +3.5 to +3.8; in the authors' CSNK2A1 substrate list)
  - SF3A1 S451
  - ZC3H8 S77
  - PAK2 T169 (weak)
- Not quantified or absent: U2SURP S37, TRAF4 S426, PHAX S368, NFXL1 S331, GPATCH2 S195, CEP152 T1241, YTHDC1 S308, CDV3 S106.

**Not obtained:**
- Borgo 2023 Open Biol, CK2α/α′ KO C2C12: HTTP 403.
- Gyenis 2022 Front Mol Biosci: zip truncated.
- Franchin: paywalled; not attempted.
- Wells 2021 (the SGC-CK2-1 probe paper) appears to contain no phosphoproteomics.

**Verdict.** The decisive test audit D asked for now exists and is positive. KPNA4 S60 is a **CK2α (CSNK2A1) catalytic-activity-dependent site**: a selective probe lowers it, and a drug-resistant allele rescues about 80% of the loss. Audit D's caveat "medium confidence it is CK2 rather than an off-target" should be raised to **high confidence**.

What remains open:
- Direct versus indirect phosphorylation. The S-x-x-D/E acidic motif makes direct CK2 phosphorylation likely, but this is not proven.
- The "fast-turnover (≤1 h)" property. The earliest time point is 4 h.
- Generality beyond U2OS.

## Files

- Scripts in `/workspace/scripts/audit_F/`:
  - `lib_f.py`
  - `01_sites.py` (a)
  - `02_akt_rsk.py` (b)
  - `03_null.py` (c)
  - `04_combo.py` (decryptM combination)
  - `05_spotcheck_menyhart.py`
  - `06_final_counts.py`
  - `t2_*.py` (Task 2 retrieval)
- Tables in `/workspace/results/audit_F/tables/`:
  - `conditions_classes.tsv`
  - `site_condition_long.tsv`
  - `site_lab_summary.tsv`
  - `site_lab_counts_final.tsv`
  - `motif_percentiles.tsv`
  - `akt_vs_mapk_verification.tsv`
  - `akt_verified_conditions_per_site.tsv`
  - `roux_rsk_conditions.tsv`
  - `decryptM_A549_combo.tsv`
  - `null_tests.tsv`
- Task 2 table: `/workspace/results/audit_F/task2_menyhart2023_module_sites.tsv`.
- Downloads (untrusted, processed tables only): `/workspace/cache/audit_F/dl/`.
