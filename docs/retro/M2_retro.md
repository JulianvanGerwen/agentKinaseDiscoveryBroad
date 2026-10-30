# M2 outside retro (fresh reviewer; read PROMPT, TREE, LESSONS, M1 retro, the 4 audits, a few SUMMARYs; no code/logs)

## 1. Stop doing
- **Stop feeding Finding A.** After the M1 retro said to stop, I22, I23 and I25 still ran (three more branches, all negative or reframed). A is now "a reproducibility asymmetry with no mechanism, effect LOR ~0.3, MEK-dependent". Report it once, as a methods caution. No more A branches.
- **Stop promoting before audit, and leaving candidates unaudited.** F (27 PSP "corrections") and G (EGFR-mut LUAD ERK ceiling) have no independent audit. F's own summary says only 7% of edges are testable and gold edges are only 22% "supported", so "unsupported" does not mean "wrong". G has no replication yet and only a PREREG.md in b_I29.
- **Stop writing tumour claims without a pre-analytical control.** The ischaemia confound (audit_I10) was found by an auditor after B was called "not ischaemia-driven". E and G were checked only afterwards (V-E).
- **Stop treating a PI-run "validation" as independent.** V-I03 said CK2 to KPNA4 S60, I26 contradicted it, and only audit_D resolved it. Validation by the producing side is where overreach starts.
- **Stop reasoning from family-level scores to named kinases.** D's calls rest on inhibitor responses, and RRAGC S15 on about 1.5 studies. LAMTOR1 "2 studies" is one lab, and the S98/S97 localisation was misattributed.
- **Process**: sub-agents still cannot write SUMMARY files (I28, Q8 and I29 have none, with ~1 h left). Use a final-message template, saved by a trivial script.

## 2. What worked well
- **The independent audit gate.** Four audits changed wording each time: the I19 null bug, I10's ischaemia confound, I24's "non-transitivity" which failed as worded, and D's S97/S98 and same-lab corrections. Each audit left a defensible sentence.
- **Pre-registered negatives**: I13, I23, I24.1 (powered, 0/162), I25, I27, I01. They are honest and cheap, and they narrow the claims.
- **The KPNA4 audit.** It resolved an internal conflict with a motif test, an independent mitotic-HeLa dataset and an off-target motif exclusion.
- **Cross-branch positive controls** (lesson 13): I26's CK2 control exposed V-I03.
- **A well-chosen null per result** (matched-null B, size-matched sets, mirror test).
- **Rebalancing to named entities** (I03, then D).

## 3. Were the M1 recommendations followed?
| Rec | Result |
|---|---|
| Auditors assign C; [!] only with audited wording | Followed for A, B, D, E, which gave C 3-3.5. Not for F or G (producer scores only). The tree still carries [!] for E and D. |
| Stop A / rebalance to named discoveries | Half. I03 and I13 launched. I03 yielded D (named sites); I13 was negative. But I22, I23 and I25 consumed A-related effort. |
| Gate I10 (adjust for subtype, proliferation, purity) | Done, but the audit found the missed confound (ischaemia). It halved B's MAPK axis. DNA damage is clean. |
| ATP "consistency check" wording | Followed (secondary). The recommended extension (more drugs, ChEMBL, decryptM 2.0) was not done, so it remains n=5. |
| Literature per claim, MCP tools | Partly. Audit_D has a literature section. There is no per-finding table, and PMC was blocked. |
| Tumour-vs-NAT, ranked top-50 edges | Q8 still running. I03 gave 4 named edges, not a top-50. |
| Cross-validation as step 1 | Mostly not. Validation still follows the branch (V-I03, V-E). |

## 4. Ranking the candidates for the final report
1. **D3 KPNA4 S60: fast-turnover CK2-type site, correcting a published MAPK attribution (C3.5/N3/P3.5).** Present as a **discovery**, with the caveat "CK2 most likely; awaits SGC-CK2-1". Best-supported and most checkable.
2. **E (narrow form): abundance-driven receptor pY does not reach ERK in tumours, while acute perturbations propagate (C3/N3/P4).** Present as a **discovery with the audited wording**. It is robust to ischaemia (V-E). The caveat is that it is "largely known" in spirit (buffering/feedback), and that only 7 cohorts are usable.
3. **D1 LAMTOR1 S97/S98 and D2 RRAGC S15: candidate mTORC1 sites (3/4/3 and 2.5/4/3).** Present as a **caveated observation** (candidate sites, not direct substrates). The evidence is one lab plus decryptM, S6K was not excluded, and AKTi sensitivity is high. Present D2 as weaker, with D1.
4. **G: EGFR-mut LUAD, lower ERK but higher MAPK transcriptional output.** Present as a **discovery only if** I29 replicates in an independent LUAD cohort. Otherwise it is a caveated observation. It has the best "disease" story, but it is unaudited and based on one cohort. The wording (ERK-toxicity ceiling) is interpretation.
5. **F: 27 PSP corrections.** Present as a **caveated observation**. Report the 3-4 literature-confirmed ones (CHEK1 S280, EPHA2 S897, CTTN S405) as validation of the method, and the new ones (DBNL S269, TFPT S180, AKAP11 T1100) as candidates. Do not claim "27 corrections".
6. **B: perturbation footprints are coherent in all 10 cohorts.** Present as a **caveated observation**. Lead with the DNA-damage footprint, which is clean. Present the MAPK axis as ischaemia-confounded. The "4x PSP" comparison is not novel as a framing.
7. **A**: report as a **methods caution** (one paragraph), not a discovery.
8. **Secondary**: ATP consistency (I17) is a **caveated observation** (n=5). EPN1 S435 (D4) is a pathway-level **observation**. **Drop**: the I24.1 offset/gain null is a one-line negative. Also drop the I10.3 immune-PTP, PAK4, Hippo, dark kinases, reporters, allele and phosphatase items (negative results only).

**Single most important remaining check.** An orthogonal check on the CK2 call and the mTORC1 sites. There is nothing in-house to settle it, so the best step is a pre-write-up literature and public-data check: does KPNA4 S60 move with selective CK2 inhibition (SGC-CK2-1, or CK2 rescue data, in PRIDE/PhosphoSitePlus)? If that is unavailable, the most important check is **G's independent-cohort replication with an ischaemia adjustment**, since G is the only candidate that could become a disease-relevant headline and the only unaudited one. Send G, with F, to an independent auditor before write-up.

## 5. LESSONS
**Merge or delete**
- Merge 1, 2 and 13 into one rule: "null next to every claim; an independent auditor recomputes it; cross-check against other branches' positive controls."
- Merge 3 into 2 (validation by an independent party first). Delete 4 (obvious). Merge 6 into a rule about data access and 12 into a rule about thread caps (OMP and similar).
- Merge 9 into 14 as "choose exposures that vary and adjust for pre-analytical variance".
- Keep 5 as an infrastructure item, but fix it rather than reminding.
- Keep 7, 8 and 10. Rule 11 is violated in practice, so make it a hard cap (at most 1 branch after a positive).

**Add**
- A claim's replication count is the number of independent labs, not studies or conditions (D1: "2 studies" = 1 lab).
- Check site localisation and identity across datasets before comparing evidence for one site (the S97/S98 mismatch).
- "Not X-driven" must be tested, not asserted, and its test script must exist (B "not ischaemia" had none).
- An inhibitor needs a verified on-target positive control in the same experiment before reading off-target or no effect (CX-4945 and CK2, S6Ki, I26's design).
- Freeze the claim list and wording at an audit date. Anything unaudited by T-1h goes in the report as an observation. Plan the audits first.
- Opposite-sign or powered nulls get a one-line entry in the tree, with no new branches.
