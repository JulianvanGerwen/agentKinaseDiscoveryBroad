# FINAL fact-check of REPORT.md and retro

Reviewer: independent, read-only. REPORT.md was not edited. Each entry reads: REPORT line (quote) → problem → correct value (source).
Method: every cited path was tested for existence. Every number was traced to `results/audit_*/AUDIT.md`, `results/val_*/SUMMARY.md`, `results/b_*/SUMMARY.md` or `RESEARCH_TREE.md`.

**Paths:** all paths in REPORT.md exist (RESEARCH_TREE.md, LESSONS.md, LITERATURE.md, docs/lit/novelty_checklist.md, every `results/...` directory, and all 7 key-figure PNGs). No wrong paths. Three citation gaps are listed under C.

## A. Mismatched numbers and scores

1. L3/L17 "Every claim below quotes the wording of an independent audit" and "scored ... by an independent auditor who recomputed it" → false for several scores.
   - K1 "4 / 3 / 3.5": audit_D gives **3.5 / 3 / 3.5**. audit_F gives no C/N/P. The "C4" is the PI's own entry in RESEARCH_TREE D9.
   - K2 "4 / 3 / 3": audit_I10 gives only the verdict "QUALIFIED" and no scores. The numbers are b_Q8's **self-scores** (preserved wiring C4/N3/P3). The same b_Q8 summary scores the DNA-damage gain **C4 / N2 / P3**, and the report gives that gain N3. b_Q8 was never audited.
   - K3 "3.5 / 4 / 3": audit_I19 says "Auditor score **C3**/N4/P3" (self-assessment C3, N3–N4). The tree I19.1 line says the same. The 3.5 comes only from tree D9.
   - K4 "3 / 3 / 4": the tree attributes these scores to the auditor, but audit_I24.md contains no numeric score. It does say the claim is "largely known" and "expected from these sources".
2. L15 "Key finding: requires C ≥ 3.5 after audit" → inconsistent with the project's own rule (tree "Ranking scheme": "key only if **C≥4**, N≥3, independently audited"). It also fails for the report's own table. K4 has C3 and K3's audited C is 3, both below 3.5. K2's "4" is unaudited. Only K1 (audited 3.5) passes the report's threshold, and none passes the tree's.
3. L16 "7 audits" → there are **6** audit directories: audit_D, F, H, I10, I19, I24. val_B and val_E are validations, not audits. The tree also says 7. "24 branches" could not be reproduced: 18 `b_*` plus 4 `val_*` is 22 (24 only if the two infra directories are counted).
4. L57 "Already present in normal tissue (Q8)" → numbers match b_Q8 (6 NAT cohorts, 91%, 529 edges). The source is an unaudited producer summary.
5. L62 "uncorrelated with the footprint scores (|ρ| ≤ 0.09)" → val_B T4 pooled ρ: MEK 0.09, EGF −0.02, TNF −0.02, ERK loop −0.05, but **DNA damage +0.139 [0.066, 0.211], p = 2e-4**. DNA damage is exactly the footprint K2 reports as the tumour-specific gain. The |ρ| ≤ 0.09 statement holds only if DNA damage is left out.
6. L58 "survives composition, ischaemia and proliferation adjustment" → b_Q8 (scripts/b_Q8/common.py, `ischaemia_fp.tsv`) adjusted for the compendium footprint **proxy**, not measured GDC times. val_B shows that proxy does not track time and over-adjusts; LESSONS 5 says to use measured time. The DNA-damage gain has not been tested against measured excision-to-freezing time, and val_B shows a weak positive correlation (see item 5).
7. L91 "Robust to ischaemia adjustment (−0.16)" → val_E: −0.158 [−0.254, −0.062]. The number is right. The covariate was the same footprint proxy; val_E did not use measured time.
8. L95 "Power was 84–94%" → b_I24 says 84–94% only for spike-ins in two lesions (DUSP6, PTEN). Not shown for the other lesions. The 0/162 result (min q .99) comes from the unaudited b_I24 and is not covered by audit_I24, which audits only candidate E.
9. L78 "decreases do not (Δ −0.70 [−1.21, −0.25])" → the Δ is dLOR for the **>2 |log2FC| bin only**. The 1–2 bin is −0.38 [−0.68, −0.16] (audit_I19 §2f). Decreases in stimulations replicate weakly (LOR 0.48 / 0.66 / 0.59, above base), so "do not" is stronger than the audit. The audit warns against the phrase "do not replicate" for the inhibitor case.
10. L110 "(split by lab, 0.27 vs 0.15 for random sites)" ✓. L111 "Fast sites are better reporters in tumours" → audit_H H3: "survives at **about half the size**" (gap .117 → .063 after ischaemia plus 5 PCs). Family-specific excess is only about .04–.05, and random non-substrate "fast" sites also out-couple.
11. L112 "+0.035 ... holds when whole labs are left out" → audit_H H4: +.035 holds only in the "edited" subset (16 studies), and that subset was not pre-declared. **All studies: +0.010 [0.000, 0.023], ns.** Condition-weighted +.015; decryptM (107 of 179 edited conditions) gives −.007. The report omits all of these qualifications.
12. L119 "did not replicate in CPTAC Confirmatory or in the Taiwanese cohort" → b_I29: Taiwan pERK +0.17 [−0.64, 0.97] is **underpowered**, not a failed replication, and Taiwan ETV4 is +0.55 (p .05), directionally consistent. Only Confirmatory (−0.03 [−0.55, +0.48]) is a clean non-replication, and that cohort has only 8 pTpY tumours.
13. L105 "Both supporting studies come from one lab" → audit_D: Hijazi and 26060313 are both Cutillas, but decryptM (Kuster) carries the dose data, at S97. With the S97/S98 merge it is "3 studies from 2 labs". The report's "matches PRAS40 S183" uses decryptM S97 data, so the line mixes the two framings. The audit's wording: "(one phosphopeptide, localisation ambiguous) ... falls in 2 labs and 4 cell lines".
14. L113 "Dasatinib engages p38α/MK2 in cells at about 1 µM" → tree I17.4 and b_I17: "weak shuffle separation", C3/**N2**/P3. It is stated without a caveat, and the ATP item has no independent audit.

## B. Claims stronger than the audited wording

- **L6/K1 L36–37 (CK2 off-targets):** "CX-4945 off-target motifs are **excluded**" → audit_D: "**depleted**" (CLK1 OR 0.5, DYRK1A 0.2, HIPK2 0.4, CDK1 0.2). It also says "No specific off-target is supported", not that all are excluded.
- **K1 L32/L35 (14-site fast set):** "All 6 measured members ... behave the same way" → audit_F: PAK2 T169 is "(weak)". The 14-site set was defined on CX-4945 specificity, and 8 of 14 sites were not quantified in Menyhart.
- **K1 omissions:** audit_D's skeptic point (Menyhart: only about 15% of CX-4945-down sites at 4 h are truly CK2-dependent). CPTAC co-variation adds nothing beyond the acidic motif (58th percentile among CK2-motif sites). audit_F also says the Menyhart arm "upgrades audit D" and does not test the ≤1 h kinetics. The report does state the 1 h limit, which is correct.
- **K1 L38 "points to CK2 inhibitors (silmitasertib is in clinical trials), not MEK inhibitors":** no audit or lit file supports the trial status or the "head-and-neck cancer" and "nuclear import" description of Hazawa 2019. audit_D says only "the only functional paper ... ('MAPK-dependent', 24 h U0126, pan-pSer pulldown)". Hazawa appears in no LITERATURE.md or checklist row. The therapeutic inference is beyond what was audited.
- **K3 L82 "The reverse layer is not a property of the cell line":** b_I22: reverse profiles are not reproducible across labs for the same cell line (MCF7 LOR −0.10), but **within Hijazi they group by cell** (r .15–.18 vs .06–.09). Correct wording: "study×culture-specific, not shown to be cell-intrinsic". audit_I19 also says it is the cell line, not the study, that opens the gap, and warns against "lost when the study changes".
- **K3 L24, L76 "portable across labs and cell lines":** audit_I19 §2e shows that cross-study replication drops for **both** directions when the cell line changes (down LOR 1.06 → 0.53). The audited statement is "decreases lose less". Replication was measured across **studies**, not labs. audit_I19 also says most portable down changes are shared across inhibitor classes (convergent nodes).
- **K3 L83 "Explanations ruled out ... All were pre-registered tests, all negative":**
  - Timing (b_I05) has no preregistration.
  - b_I22 used "pre-declared" H1–H4.
  - I25 lists deviations from PREREG (added unmatched control, T3d, T2b).
  - I23 was **refuted in sign** and is "fragile ... underpowered" (about 8 informative drugs, 3 lines).
  - b_I25 states that "shared-control / within-run structure" is untested, and I22 states H3 "saturation untestable".
  - "Ruled out" is too strong; "not supported" is the summaries' wording.
- **K3 L85 "Compensatory or feedback increases ... fail to transfer":** audit_I19 lists "reverse = compensatory" as "interpretation, not shown". The inhibitor arm is MEK-dependent: leave-MEK-out ΔLOR 0.20 [−0.05, 0.53] includes 0, and the effect is absent at ≤30 min. DDR inhibitors show the reverse. The report mentions only "leans on MEK".
- **K3 L80 "In tumours: only PI3K-axis down-sites predict PIK3CA/PTEN (AUROC 0.56–0.57)":** audit_I10 gives .564 (p = .01), size-matched .554, up .460, "real but small". audit_I10 also says the classifier result "does not replicate Finding A as a general down-site property" but supports primary direction. The report calls it "independent replication" in the K3 table.
- **K2 L51–54 "This survives ... removal of the top 5 PCs":** audit_I10: PC removal cuts MAPK/EGF/TNF coherence by about 55–65% (MEKi .153 → .068). Only DNA damage and metformin are unaffected. The audit's own wording is "remain coherent". The halving is not reported.
- **K2 L59–63:** the audit's recommended wording says the MAPK axis "tracks ERK activity / MAPK1/3 ρ≈0.8" and "adds nothing over kinase scores for RAS/RAF". The report mentions "do not beat kinase scores for genotype" only under Limits.
- **K4 L91 "Receptor abundance lowers ERK":** audit_I24 says "goes with" (a cohort-level association, not a causal effect). The Fritsche-Guenther link is the auditor's plausibility argument.
- **K4 L94 "Measurement noise is excluded":** audit_I24 says classical noise is rejected by the simulation (P = .010). Correlated (non-classical) error is "the remaining measurement route", partly proxied by global pY, and the simulation "assumes independent errors". Correct wording: "classical measurement noise does not explain it".
- **K4 audit caveats omitted:**
  - The −0.15 comes from the unadjusted `tc` variant. On protein-adjusted data the gap is not significant (p = .16).
  - The result depends on the ERK isoform (ERK2-only p = .10).
  - "10 cohorts" is really 7.
  - Only LSCC is individually significant in the IV decomposition.
- **L122 "every key finding was checked ... by logged web searches":** there is no per-finding literature table. audit_D "could not open full texts (PMC blocked)". The "logged" searches were not found. The M2 retro already noted that no per-finding table exists.
- **L143 (§4):** "PAK4 S148" is listed among negatives. The M1 retro and tree V1 say "unreplicated/untestable", not refuted. This is acceptable only because the heading says "abandoned".

## C. Citation errors and gaps

- K3 Files (L87) omit audit_I10 and b_I10 (the classifier and tumour PI3K claims) and b_I05/b_I22/b_I25 are cited, but the "3 audits" in the K3 table cannot be traced. Only audit_I19 audits the K3 claim, plus the audit_I10 §B4 direction section.
- K2 Files (L72) omit val_E, which holds the ischaemia-robustness result for K4. K4 Files (L100) list val_E but not val_B, which supports L62.
- K2 L50 "median AUROC 0.71; EGF 0.96, DNA damage 0.97, MEK 0.88" comes from b_I10 (unaudited). audit_I10's own LOSO gives inhibitors down .65 vs up .55, stimuli up .78 vs down .55 (28 footprints).
- L113 and L112 give audited numbers for audit_H and audit_F correctly: F "2 / 2 / 2.5" ✓ (audit_F: C2/N2/P2.5), H "3 / 2.5 / 2.5" ✓. The report omits the audit_F finding that "PSP already lists RSK for EPHA2 S897 and CHEK1 S280, and ERK for CTTN S405" (L116 lists only DBNL and AKAP11 as new, which is fine) and that the CDK→ERK part is "at chance".
- L103–L119: the b_I08 numbers (81/117) are unaudited producer numbers, and I08.3 reporter sites were rejected by V1 (the report does not claim them, fine).

**Verified correct (spot list):** KPNA4 z-values and rank 63/12,455; 0/25 ERK units; OR 16, p = 1.5e-6; 64%; SGC-CK2-1 −3.73 (FDR 3e-4) and +2.98 (FDR 4e-4), about 80%; Hazawa 24 h U0126; GDC n = 723 in 7 cohorts, 3 cohorts without times; coherence about 4× PSP; z ≈ 20–30; Q8 numbers; −0.15 [−0.25, −0.05], p = .002, 6/7; 0.31 vs 0.44, p = .16; 0.58 vs 0.51; P = .01; 0/162; KRAS +0.6; the LAMTOR1 0.11th percentile; RRAGC 6/6; the PSP audit 7% / 17% vs 5%; ATP 1.5–3.6 mM, n = 5; the audit_H lab split .27 vs .15; the 14-site count; the Hazawa and Menyhart bibliographic details.

---

# Final retro (≤½ page)

**(a) Stop doing.**
- Stop letting producers or the PI score C/N/P. Three of four key-finding scores in the report are not auditor numbers, and the report's key-finding bar (C≥3.5) differs from the tree's (C≥4).
- Stop adjusting for signature proxies. The "ischaemia footprint" was used again in Q8 and val_E after val_B showed it is not a time measure, and the DNA-damage gain still carries a +0.14 correlation with measured time.
- Stop running audits that cover only part of a result. b_Q8, b_I24.1, b_I17 and b_I08 were never audited, yet were written up as findings.
- Stop writing the report from tree summaries. Write it from the audit text and have each sentence checked against it. The errors above are mostly omitted caveats and "ruled out" language.
- Stop feeding Finding A. Three more branches ran after M1 said stop.
- Stop the SUMMARY-write workaround. It cost PI context twice.

**(b) What worked.**
- The independent audit gate changed the wording of every key finding, and the audits stand as the report's reliable core.
- The audit_D to audit_F chain, with a public chemical-genetic dataset fetched within the audit's time box, turned a conflict into the one strong result.
- Pre-registered negatives (I13, I23, I24.1, I25, I27) were cheap and narrowed the claims.
- Measured GDC times (val_B) corrected an audit's over-reading.
- Cross-branch positive controls exposed V-I03. Concurrent branches stayed within the quarter-node cap. Honest dropping of G.

**(c) Does the outcome meet the brief?** Partly.
- *Discoveries:* one reasonably compelling, named-entity discovery (K1, audited C3.5/N3). K2–K4 are audited at C3, "largely known in spirit", or methodological. "Several sufficiently novel and compelling" is therefore not clearly met.
- *Diversity:* met. About 24 branches across the compendium, CPTAC, dose pharmacology, causal IV, phosphatases, dark kinases, turnover and pharmacology. Positives cluster in compendium meta-statistics and CPTAC coherence.
- *Literature corroboration:* partial. There is a 97-row checklist and good audit-level checks (Menyhart, Banks, Hazawa), but no per-finding table, PMC was blocked, and the therapeutic claim in K1 has no source.
- *Ranking:* done but inconsistent. Scores are mixed producer and auditor, and the thresholds disagree.

**Most valuable next step.** Test K1 directly, by a human or in a wet lab. This means an in vitro CK2α phosphorylation assay on KPNA4 S60, plus a ≤1 h SGC-CK2-1 time course with the inhibitor-resistant CK2α rescue, ideally with nuclear-import read-out. It is the only finding where direct kinase attribution, the ≤1 h turnover claim and the Hazawa correction are still open, and computation has saturated. In parallel, apply the corrections above to REPORT.md and have a human review it, since "no result has been verified by a human".
