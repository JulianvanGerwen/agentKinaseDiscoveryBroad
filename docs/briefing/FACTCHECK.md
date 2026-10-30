# Fact-check of docs/briefing/BRIEFING.md

Reviewer: independent, read-only (BRIEFING.md and figures not edited). 30 Oct 2026.
Each item: **location (quote / line)** → problem → suggested fix. Line numbers refer to BRIEFING.md as checked.
Sources were used in this order of authority: `results/audit_*/AUDIT.md` > `results/val_*/SUMMARY.md` > `results/b_*/SUMMARY.md` > `RESEARCH_TREE.md` > `REPORT.md` / `docs/retro/FINAL_factcheck_retro.md` > `LITERATURE.md` / `docs/lit/novelty_checklist.md`. Citations were checked against PubMed.

Overall: most numbers are right and many corrections from the earlier fact-check are respected (score attributions, 6 audits not 7, PC-halving, ischaemia proxy vs measured time, H4 subset, ERK-ceiling replication wording). The problems below are mostly (i) process claims that overstate audit coverage, (ii) a few earlier fact-check corrections that were not carried over, (iii) statements the figures themselves contradict, and (iv) five reference-list errors.

---

## CRITICAL

**C1. Audit coverage is overstated (process claims).**
- L15 "It also sent every candidate claim to an independent auditor" and L38 "Independent audits changed the wording of every positive claim and killed two"; L341 "Every positive claim was narrowed by its audit, two candidates were killed".
- Problems:
  - Several claims in the briefing were never audited: the tumour-vs-normal comparison and DNA-damage gain (b_Q8), "offset not gain" (b_I24 finding 3; audit_I24 audits only candidate E), the ATP/Cheng–Prusoff check (b_I17), and activity vs abundance (b_I08). The earlier fact-check retro says this explicitly ("Stop running audits that cover only part of a result. b_Q8, b_I24.1, b_I17 and b_I08 were never audited"). The briefing itself admits Q8 was unaudited (L213), so L15 contradicts L213.
  - Of the "two killed", only immune–phosphatase mimicry failed an audit (audit_I10 verdict C). The ERK ceiling (G) was dropped because its own pre-registered replication failed in b_I29, not because of an audit. RESEARCH_TREE says "dropped after replication or audit".
  - K1 was not simply "narrowed": Fig 8 shows it upheld and then strengthened. "Every positive claim was narrowed" contradicts Fig 8.
- Fix:
  - L15: "sent the main candidate claims to an independent auditor".
  - L38/L341: "Independent audits changed the wording of every key finding they covered. One candidate failed its audit (immune–phosphatase mimicry) and one failed pre-registered replication (ERK ceiling). Some components (tumour vs normal tissue, offset-not-gain, ATP check) were never audited."

**C2. Measured ischaemia time vs footprint scores (earlier fact-check item 5 not respected).**
- L206 "measured time was unrelated to the footprint scores (pooled |ρ| ≤ 0.09)".
- Problem: val_B T2/T4 gives the DNA-damage score a pooled ρ of **+0.14 [0.07, 0.21]**. |ρ| ≤ 0.09 holds only for MEK, EGF, TNF and the ERK loop. The briefing's own caveat (L214, "ρ = +0.14") contradicts L206, and DNA damage is the footprint behind the tumour-specific gain.
- Fix: "measured time was unrelated to the MEK, EGF and TNF footprint scores (pooled |ρ| ≤ 0.09); the DNA-damage score showed a weak positive correlation (ρ = +0.14)".

**C3. "every core link has a positive gain in all 10 tumour cohorts" (L253).**
- Problem: audit_I24 §0: "'10 cohorts' is really 7". The receptor (RTK) stage is missing in OV, PDAC and UCEC, and in BRCA and COAD it is a single site. The briefing's own Fig 6b shows (7) for receptor pY → adaptor pY and for mTORC1 → S6.
- Fix: "every measured core link has a positive gain in every cohort where it can be measured (7–10 cohorts; Fig 6b)".

**C4. "Offset, not gain" is stated more strongly than the evidence (earlier fact-check item 8 not respected).**
- L33 "Cancer mutations … do not change the pathway's sensitivity"; L259 "None changed a gain … so this is an informative negative".
- Problems:
  - b_I24 shows 84–94% power only in spike-ins for two continuously measured regulators (DUSP6, PTEN).
  - For mutations the intervals are wide. In Fig 6c the KRAS gain change is about +0.25 with a CI reaching about +0.6, and PTEN mutation about +0.18; both point estimates lie outside the ±0.15 "detectable" band.
  - The result comes from the unaudited b_I24 (audit_I24 covers only candidate E).
- Fix:
  - L33: "no lesion detectably changed the pathway's sensitivity (well powered for regulators such as DUSP6 and PTEN, not for individual mutations)".
  - L259: limit "informative negative" to the continuously measured regulators.
  - Add to the §5.3 caveats: "the offset/gain analysis was not independently audited".

**C5. Reference-list errors (verified on PubMed).**
- **Abelin 2016**: the journal is *Mol Cell Proteomics* 15:1622–41, not *Cell Syst* (PMID 26912667 is correct).
- **Invergo 2020**: PMID 32437683 is "Prediction of signed protein kinase regulatory circuits", *Cell Syst* 10:384. "Sub-minute phosphoregulation" is a different paper (Plasmodium). Delete that phrase.
- **Hazawa**: the entry's title is invented. The real reference is Hazawa M et al., "Disease-specific alteration of karyopherin-α subtype establishes feed-forward oncogenic signaling in head and neck squamous cell carcinoma", *Oncogene* 39:2212–2223 (2020; online Dec 2019), PMID 31822798, PMC7056645.
- **Koenig 2026**: the authors are Koenig C, Cho H, Emdal KB, Piga I, Sabatier P, Lozano-Juárez S, Martinez-Val A, Olsen JV (cite as "Koenig C et al."). Title: "A phosphoproteome atlas of human cell lines reveals the landscape of kinase activity", *Nat Struct Mol Biol* 33:1561 (PMID 42816593 is correct). LITERATURE.md line 93 has the same author error.
- **Müller-Dott 2025**: the journal is missing. It is *Nat Commun* (2025), "Comprehensive evaluation of phosphoproteomic-based kinase activity inference", PMC12098709.

**C6. Fig 8 implies a K1 score that no audit gave.**
- Fig 8, K1 row: a dashed arrow labelled "then raised further" extends from 3.5 to about 4.3.
- Problem: audit_D gives C 3.5; audit_F raised *confidence* but assigned no C/N/P. The earlier fact-check (A1) flagged the "C4" for K1 as the PI's own number. Placing the arrow on the score axis presents an unaudited score of about 4+. It is also inconsistent with L183 ("audited score is Compelling 3.5").
- Fix: end the K1 arrow at 3.5 and show audit_F as a text annotation ("confidence raised to high by public rescue data; no new score"). Add to the briefing caption that non-numeric verdicts are mapped to numbers (the figure footnote says so, the caption does not).

---

## IMPORTANT

**I1. mTORC1 candidates: earlier fact-check item 13 not respected, and the Fig 7 caption contradicts the plot.**
- L287 "with dose potency matching a known direct mTORC1 site"; L289 "both LAMTOR1 studies come from one laboratory"; Fig 7a caption "LAMTOR1 S98 falls to the 0.11th percentile … close to the direct mTORC1 site PRAS40 S183".
- Problems:
  - audit_D: the dose data (decryptM, Kuster lab) exist only for **S97**. With the S97/S98 merge the evidence spans "3 studies from 2 labs". The two Cutillas studies give one lab for the S98 call only.
  - The "potency equals PRAS40 S183" comparison is the MK2206 (AKT inhibitor) fit, which is "not significant", and audit_D says it does not show direct phosphorylation.
  - In Fig 7a, LAMTOR1 S98's median z under catalytic mTOR inhibitors is about −5.0 and PRAS40 S183's is about −10.6, so "close to" is visibly wrong.
  - Text omits that LAMTOR1 S98 is unusually AKT-inhibitor-sensitive (AKTi/mTORi ratio 0.35, vs 0.09 for PRAS40) and that RRAGC S15 rests on about 1.5 independent studies (audit_D D2).
- Fix:
  - Use audit_D's wording: "LAMTOR1 S97/S98 (one phosphopeptide, localisation ambiguous) … falls in 2 labs and 4 cell lines; S6K, SGK and other downstream kinases are not excluded".
  - Drop "dose potency matching" or say "a non-significant AKT-inhibitor dose fit has potency similar to PRAS40 S183".
  - In the caption, replace "close to PRAS40 S183" with what the plot shows.

**I2. Finding 5.2: "laboratories" should be "studies".**
- L32 "reproduce across laboratories"; §5.2 heading "travel between laboratories"; L247 "the cross-laboratory test".
- Problem: audit_I19 measured replication across **studies** (41 inhibitor and 52 stimulation studies), many from the same lab (Cutillas, Kuster). The earlier fact-check (B, K3) made the same correction. audit_I19 §2e also finds that changing the **cell line**, not the study, opens the gap. The only explicit cross-lab test is one MCF7 comparison (b_I22). Fig 5 correctly says "studies".
- Fix: say "across studies" throughout, and mention the single MCF7 cross-lab comparison separately.

**I3. Finding 5.2: audit caveats omitted, and "ruled out" retained.**
- L237 "Several pre-registered or pre-declared tests ruled out obvious explanations"; L245 "leans on MEK inhibitors".
- Problems:
  - The earlier fact-check said "ruled out" is too strong; the summaries say "not supported". Timing (b_I05) was not pre-registered. I23 was "refuted in sign", "fragile" and underpowered. I25 lists deviations from its pre-registration, and its shared within-run structure is untested.
  - audit_I19 caveats not mentioned:
    - leaving out MEK roughly halves the inhibitor effect, and the CI includes 0 (0.20 [−0.05, 0.53]);
    - there is no asymmetry at ≤30 min;
    - DDR inhibitors show the reverse;
    - most portable decreases are shared across inhibitor classes (convergent nodes).
  - L235 "only the primary-direction sites of the PI3K footprint predict PIK3CA/PTEN": audit_I10 gives AUROC 0.56, "real but small".
- Fix:
  - Change "ruled out" to "did not support".
  - Add one caveat sentence covering MEK dependence (CI includes 0 without MEK), no asymmetry at ≤30 min, and the DDR reversal.
  - Add "(weakly; AUROC 0.56)" to L235.

**I4. K1: fast CK2 set and Fig 3 caption.**
- L176 "All six measured members of our fast-responding CK2 group behaved the same way"; Fig 3 caption "The five other quantified members … behave the same way".
- Problems:
  - audit_F lists PAK2 T169 as "(weak)". 8 of the 14 sites were not quantified.
  - In Fig 3b, under CX-4945 at 4 h one module site is *not* rescued (it falls further in the resistant line). Under SGC-CK2-1 at 4 h only two other module sites are plotted.
- Fix: "all six quantified members fall under SGC-CK2-1 and are at least partly rescued (one weakly); 8 of 14 were not measured".

**I5. K1: the MAPK correction is wider than the evidence, and Hazawa is presented incompletely.**
- L29 "never under ERK inhibition"; L181 "The 24-hour MEK effect in the original paper is probably indirect".
- Problems:
  - audit_D: the 0/25 result "addresses only acute (≤2 h) dependence; Hazawa used 24 h U0126". "Probably indirect" appears in no audit.
  - Hazawa et al. also used S60A/S60D mutants (functional effects on RREB1 binding and growth), and the site was first reported by Olsen 2006. L141 calls it "the only functional study" but describes only the U0126 experiment.
  - "Serine 60 lies in the region that binds importin-β" (L141) has no source.
- Fix:
  - L29: "not under acute (≤2 h) ERK-pathway inhibition (0 of 25 experiments)".
  - L181: "the 24-h MEK effect was not tested here and may be indirect".
  - Mention the S60A/D data so the critique is fair.
  - Source or soften the importin-β-binding statement.
  - Also say that the rescue data are at 4 h, so the "fast-turnover (≤1 h)" property of the CK2 set rests on compendium data only (audit_F "still open").

**I6. K4: audit_I24 caveats omitted.**
- §5.3 caveats (L267–273).
- Problems: not mentioned are:
  - the effect varies by cohort: LUAD carries the largest effect (−0.38), BRCA is positive, and only LSCC is individually significant in the IV decomposition;
  - the BRCA instrument is entirely ERBB2 amplification;
  - the IV gap depends on the ERK isoform (ERK1-only p = 0.001, ERK2-only p = 0.10);
  - correlated (non-classical) measurement error is "the remaining measurement route";
  - audit_I24 calls the result "largely known … expected from these sources".
- Heading L249 "receptor overexpression is buffered" is an interpretation; the audit's wording is "goes with".
- Fix: add one caveat bullet with these points, and reword the heading to "…receptor abundance does not pass to ERK…".

**I7. Turnover atlas: the tumour-reporter claim is stated without the audit's size qualifier.**
- L294 "Fast-responding sites are better reporters of kinase activity in tumours."
- Problem: audit_H H3: "survives at about half the size" (Δρ 0.117 → 0.063 after ischaemia + 5 PCs). The family-specific excess is only about 0.05, and random non-substrate "fast" sites also out-couple.
- Fix: "…modestly better reporters (Δρ ≈ 0.06 after adjustment, about 0.05 family-specific)".

**I8. Fig 4a: normal-vs-tumour comparison on PC-adjusted data.**
- Fig 4a ("protein- and 3-PC-adjusted") and its caption.
- Problem: b_Q8's "critical control" is that global PCs absorb 68–73% of MAPK-footprint variance in tumours but 34–48% in normal tissue, so "PC adjustment biases T−NAT comparisons" and "apparent RTK attenuation = PC artefact". Fig 4a shows MEK/EGF coherence falling in tumours, and a reader will read that as biology.
- Fix: add to the caption: "PC adjustment removes more variance in tumours, so the apparent tumour decrease for MEK/EGF is not interpretable; the DNA-damage gain also holds without PC adjustment (meta z 5.3)".

**I9. Inconsistent finding labels.**
- §2 table "Finding 2", "Finding 4"; L341 "the chemical-genetic rescue for Finding 4"; "Finding 5.2" elsewhere; Figs 2 and 8 use K1–K4, which the text never defines.
- Problem: "Finding 4" means K4 (cascade, I24) in the §2 table but KPNA4 (§4) in §8.
- Fix: define K1–K4 once in the Summary and use them everywhere: K1 = KPNA4, K2 = tissue programmes, K3 = portability, K4 = cascade.

**I10. Fig 7c caption: "selumetinib … shows no such shift".**
- Problem: the plot shows a shift of about +0.75 in the opposite direction (cells *more* sensitive than the ATP-free affinity). b_I17 calls this "allosteric MEKi overshoot (+0.4 to +0.8)".
- Further gaps in the text:
  - the ATP check is unaudited and rests on 3 Km values and 5 drugs;
  - MK2206 is an exception.
- Fix:
  - Caption: "…shows no ATP-competition shift; instead its cellular potency exceeds its affinity (+0.75 log)".
  - Add "not independently audited" to the text.

**I11. §7 and the negatives table overstate pre-registration and results.**
- L312 "Each was stopped on a criterion set before seeing the result".
- Problem: only I13 (positive control), I22 (pre-declared), I23, I25, I27 and I29 had pre-declared criteria. I01, I05, I06 and I20 had none (the earlier fact-check noted timing had no pre-registration).
- I27 row "Predicted signs mostly right": b_I27 shows P1 in the right direction but not significant, P2 in the *opposite* direction (and "G12C ≈ G12V" contradicted), P3/P3b unanswerable, and P4/P5 not confirmed.
- Fix:
  - L312: "Most were stopped on a pre-declared criterion…".
  - I27: "Predictions mixed (one right sign, one opposite), and CPTAC has too few tumours per allele".

**I12. Resource and thread-cap claims understate the breaches.**
- L345 "Four agents at four threads each kept the 16-thread allowance busy throughout"; L354 "One agent briefly exceeded the thread cap".
- Problem: LESSONS 11 says "I23 ran ~27 cores briefly", above the whole team's 16-core cap, not just one agent's 4. audit_H separately reports running 8 threads. "Busy throughout" has no source.
- Fix: "Two agents briefly exceeded their thread allowance, one of them (I23, about 27 cores) exceeding the team's 16-core cap", and drop "throughout".

**I13. Tier 2 definition vs what was audited.**
- §3 item 6, "Tier 2: audited and supported with qualifications … and replicated".
- Problem: K2's normal-tissue and DNA-damage parts (Q8) and K4's offset-not-gain half are unaudited, and K4 has no numeric audit score. This does not match the definition.
- Fix: in the Summary Tier-2 bullets, mark the unaudited components, or narrow the Tier-2 definition to "core claim audited".

---

## MINOR

1. **L163 order:** "In HL60 and NTERA2 cells, CX-4945 is the first and third strongest of about 58 inhibitors" → audit_D: #1 of 61 in **NTERA2** and #3 of 55 in **HL60**. Swap the order.
2. **L171** "inhibitors of those off-target kinases do not lower KPNA4 S60": audit_D notes one unexplained outlier (GSK3i LY2090314 in HL60, z −6.17, not reproduced). Add "(one unreproduced outlier)".
3. **Fig 3 caption** "falls sharply … in all three cell lines": HL60 is z −2.5 (5th percentile). Say "falls in all three, most strongly in mitotic HeLa".
4. **L177** "fell just outside their cut-off": the residual H/L was −0.82 against a ±0.585 window. Drop "just".
5. **L153** "raised the chance … about 2.5-fold": these are odds ratios (2.3 and 2.8), and b_I03 says "family-level only". Say "raised the odds about 2.5-fold (family-level only)".
6. **L194** footprint AUROCs (0.71; EGF 0.96 …) are unaudited b_I10 numbers; audit_I10's own leave-one-study-out gives lower values (inhibitors 0.65 down / 0.55 up). Label them as producer numbers, or use the audit values.
7. **L199** "The one robust tumour-specific change is DNA-damage": b_Q8 also reports a gained PSP CDK2 set (z 4.1–5.4). Say "the main".
8. **§5.1 caveats:** add val_B's caveat that recorded times span a narrow window (IQR about 8–35 min), so effects acting only beyond about 1 h are untested.
9. **L259 / Fig 6c** KRAS offset: LSCC (5 mutants) contributes about 0, so the effect is carried by LUAD and COAD (val_E).
10. **L342** "about seven minutes": the source is LESSONS 3 ("V1: 3/4 failed in ~7 min"). Fine; note that V1 also *validated* IKBKB S697, so the §7 V1 row "new kinase reporter sites … failed" should read "most reporter sites failed".
11. **L350** "the PI still ran three more branches on it, all negative": per the RESEARCH_TREE M1 status, I22 and I23 were already running at M1, and I22 produced the reframing (marked [x]). Say "three more branches (two already running) went to it; none found a mechanism".
12. **Fig 8:** G is plotted with self-score 3, but b_I29's self-score is C2. "First proposed (self-assessed)" is therefore inaccurate for G; use 2 or note the mapping.
13. **Fig 1a** "Fetched during the project" omits the Kuster KRAS dose dataset listed in L75. **Fig 2** omits I06 and I18/I21, which the text cites (§2 table, §7 "Others").
14. **L62** Kinobeads "~240 inhibitors": Klaeger 2017 profiled 243 clinical drugs.
15. **Credit gaps:**
    - The Cheng–Prusoff relation is used but uncited (Cheng & Prusoff 1973, *Biochem Pharmacol* 22:3099).
    - "Control engineering: gain vs offset" in the §2 table has no reference.
    - "CoPheeMap 2025" (L217) is missing from the reference list: Jiang W et al. 2025 *Nat Commun* 16:2766, PMID 40113755 (already in LITERATURE.md).
16. **Reference details:**
    - Yang G is *EMBO J* 38:e99548 (2019; online Dec 2018), PMID 30552228; give the PMID.
    - Bayer 2025's real title is "Chemical proteomics decrypts the kinases that shape the dynamic human phosphoproteome" ("decryptM 2.0" is a label).
    - Kanshin 2015's title is "A cell-signaling network temporally resolves specific versus promiscuous phosphorylation".
    - Kanshin 2017: PMID 28265048.
    - Unni 2015 / 2018: PMIDs 26047463 / 30475204.
    - Add the missing PMIDs:

      | Reference | PMID |
      |---|---|
      | Cundell 2016 | 27551054 |
      | Davis 2011 | 22037378 |
      | Fritsche-Guenther 2011 | 21613978 |
      | Huang & Ferrell 1996 | 8816754 |
      | Ling 2020 | 31993556 |
      | Mertins 2014 | 24719451 |
      | Arshad 2019 | 31227600 |
      | Menyhart 2023 | 41883516 |

    - Arshad, Chen 2020, Cundell, Davis, Fritsche-Guenther, Hazawa, Ling, Menyhart, Mertins 2014, Patricelli, Unni, Vasta and Yang are not in LITERATURE.md or the novelty checklist. They were never vetted at the literature stage; they are now checked here.
17. **Jargon not explained for a broad audience:**
    - FDR and q-value (Fig 3, Fig 6);
    - "PC" / "3-PC-adjusted" / "global patterns of variation" (Fig 4; L195 uses the plain wording, the captions do not);
    - "meta z", "standardised slope", "SD" units;
    - pEC50 and pKd (Fig 7c axis);
    - "spike-in simulations", "cluster-bootstrap CI", "leave one study out";
    - "activation loop", "mitotically arrested", "chemical-genetic rescue" (explained only implicitly in L173);
    - "instrumental variables" and "co-essentiality" (§2 table).

    A short glossary box would fix this.
18. **L183 / §9** "fast-turnover CK2 sites": define once in plain words ("sites whose phosphate is lost within an hour of inhibiting the kinase").

---

## Verified correct (spot list)

- **Project counts:** 1,937 / 178 studies / about 190k sites / 850 inhibitor conditions; CPTAC 1,030 tumours and 523 NAT; 18 branches, 4 validations, 6 audits, 2 milestone retros; 9 negative branches; 5 caveated observations; 18 h; quarter node (16 cores, about 94 GB, which implies a 64-core / 375 GB / 4-GPU node); decryptM 2.0 restricted, 133 drugs; about 15 overstatements and 0 fabricated numbers (LESSONS 15).
- **K1:** z −7.0, rank 63 of 12,455; 0 of 25 below −2, median +0.19; OR 16, p 1.5e-6, 64%; SGC-CK2-1 −3.73 (about 13-fold), rescue +2.98 (about 80%, FDR 4e-4); Fig 3b FDRs (7e-2, 3e-2, 2e-1); HEDICED**S**DIDGDYR; score 3.5/3/3.5; Menyhart, PXD038050, U2OS.
- **K2:** 4× PSP; PC halving; 529 relationships, 91% / 91%; meta z 6.3; n = 723 in 7 cohorts; 0.170 → 0.169 → 0.087; ρ −0.03; DNA-damage ρ +0.14; 3 cohorts without times.
- **K3:** 41 / 52 studies; LORs 0.32 → 0.77, 0.23 → 0.35, 0.56 → 1.62, 0.48 → 0.59; 12 of 16 and 9 of 12 footprints; score 3/4/3.
- **K4:** 0.10 [0.02, 0.18] vs 0.17 expected; −0.15 [−0.25, −0.05], p 0.002, 6 of 7 cohorts, n = 613, −0.16 with proxy; KRAS +0.67 (q 0.005) / +0.60; 162 tests, smallest q 0.99; 84–94% power.
- **Caveated observations:** LAMTOR1 0.11th percentile; RRAGC 6/6 organoid lines; turnover 0.27 / 0.15 / 0.17, CDK 0.10; ATP 1.5–3.6 mM, median 2.3, 1.9–2.8 log units, n = 5; activity vs abundance 81 of 117; PSP audit 7%.
- **Negatives:** dark kinases about 120 and 66%; PTEN control 10/10.
- **Citation facts:** P100 has 96 sites; Ochoa 2016 covers 215 kinases; Needham "<5%"; Ling 2020 reports AMPKα1 S347; Cundell 2016 shows the B55 pThr preference; Hazawa's 24 h U0126 and pan-pSer read-out are described accurately.
- **Paths in §10:** README.md, NOTICE.md, docs/PATHS.md, PROMPT.md, LESSONS.md, scripts/ and docs/lit/ all exist.
