# Audit D: I03 / V-I03 kinase→site nominations (LAMTOR1 S98, RRAGC S15, KPNA4 S60, EPN1 S435)
Auditor: independent re-analysis, 2026-10-29. I did not modify any I03 or V-I03 file.
Code: `scripts/audit_D/` (common.py, 01–09). Tables: `results/audit_D/tables/`. Logs: `cache/audit_D/0*.log`.

**Method.** I used a robust per-condition z (median/MAD over all sites) on the full compendium.
- Units are study × agent × cell. Doses below 100 nM are dropped, dose series keep only 1 µM, and reversed or background contrasts are excluded.
- Leave-one-study-out (LOSO) is applied to every class median.
- A genome-wide percentile is computed for each class median, against sites measured in at least as many units.
- Negative controls are the same protein's other sites and motif-matched random sites.
- Motif-set tests use Cantley percentiles (b_I05 `motif_pct.parquet`) and Klaeger Kinobeads (Table S2, Silmitasertib).
- Dose data come from decryptM `curves_all`.
- CPTAC co-variation uses protein-adjusted tumour residuals, Stouffer-combined, with the candidate's percentile taken against all sites and against motif-matched sites.

## Summary of verdicts
| Claim | Verdict | C / N / P |
|---|---|---|
| D1 LAMTOR1 S98 is a direct mTORC1 site | **Qualified**: it is a strongly mTOR-kinase-dependent site, but direct phosphorylation is not shown and the residue is ambiguous between S97 and S98 | 3 / 4 / 3 |
| D2 RRAGC S15 is a direct mTORC1 site | **Qualified (weaker)**: mTOR-dependent and the motif fits, but the evidence rests on about 1.5 studies | 2.5 / 4 / 3 |
| D3 KPNA4 S60 depends on CK2, not MAPK | **Survives with qualification**: it is a fast-responding CK2-type site, not an identifiable off-target; there is no second CK2 chemotype | 3.5 / 3 / 3.5 |
| D4 ERK pathway → EPN1 S435 | **Survives at pathway level**; direct ERK phosphorylation is unproven | 3 / 3 / 3 |

---
## D1. LAMTOR1 S98 (RLAVLSS[S]LTHWKKL)
**Reproduced.**
- **Catalytic mTOR inhibitors:** median z −5.07, 8/10 units below −2, 0.11th percentile of 5,459 comparable sites.
  - LOSO range −4.9 to −6.1.
  - 2 studies: Hijazi 2020 (31959955; HL60/MCF7/NTERA2) and 26060313 (MCF7).
  - Both studies come from the same lab (Cutillas). The "2/2 studies" replication is therefore within one lab.
  - PI-103 (dual PI3K/mTOR inhibitor) gives −7.07.
- **Same-protein negative controls are flat:** LAMTOR1 S27 +0.05 (15 units), S56 −0.66, S26 −0.07.

**Problems the V-I03 summary does not state clearly:**
1. **Localisation.** In decryptM, S98 has **no** AZD8055, dactolisib or pictilisib curve. The quoted values belong to **LAMTOR1 S97**:
   - AZD8055 pEC50 7.50 (Δ −2.08*)
   - dactolisib 7.80 (Δ −2.27*)
   - pictilisib 5.45 (Δ −2.43*)

   decryptM localises the phosphate to S97, while the Cutillas datasets localise it to S98 (S96, S97 and S98 are consecutive serines). The independent-lab dose evidence therefore applies only if S97 and S98 are the same site. The claim should read **"LAMTOR1 S97/S98 (ambiguous)"**. With that merge it spans 3 studies from 2 labs.
   - The Cantley MTOR percentile differs by residue: S98 75, S97 36, S96 21. The motif fit is weak in any case.
2. **The S6K-inhibitor negative is uninformative.** S6K inhibition is not verified in those conditions:
   - The PSP RPS6KB1 substrate set has median z −0.1 to −1.0 across the 7 S6Ki conditions.
   - RPS6 S236 rises in 4 of them (+1.3 to +2.3). RPS6 S240/S244 were not measured.
   - Only EIF4B S422 falls (−2.1 and −1.2).

   "S6Ki ratio ≈ 0" therefore cannot exclude S6K.
3. **Upstream sensitivity is above that of the direct controls.** Class-median ratios:

   | Site | AKTi/cat-mTORi | PI3Ki/cat-mTORi |
   |---|---|---|
   | LAMTOR1 S98 | **0.35** | 0.22 |
   | PRAS40 S183 | 0.09 | −0.07 |
   | 4EBP1 S65 | 0.12 | −0.02 |
   | LARP1 S766 | 0.04 | 0.11 |
   | EIF4B S422 (S6K/RSK) | 0.40 | 0.25 |
   | GSK3B S9 | 0.53 | 0.31 |

   - LAMTOR1 S98 sits between the direct mTORC1 class and the AGC class.
   - In MCF7, AKT inhibitors lower it by −3.2 to −6.3.
   - In decryptM MK2206 (A549) it has a pEC50 of 6.5–6.85, but the fit is not significant. That potency equals PRAS40 S183 (6.4–6.75). It is compatible with mTORC1 lying downstream of AKT, but does not show that mTORC1 phosphorylates the site directly.
4. **CPTAC gives no support.** Co-variation with 4EBP1 T37/T46 has meta-z −0.9 (23rd percentile, 3 cohorts).

**Skeptic's view.** All that is shown is "falls fast under ATP-competitive mTOR inhibitors". That holds equally for any S6K, SGK or other mTORC1/2-dependent kinase whose inhibition was not verified here. A residue-ambiguous peptide, one lab for the S98 call, and no rapalog data leave "direct" as a hypothesis.

**Defensible wording.** "LAMTOR1 S97/S98 (one phosphopeptide, localisation ambiguous) is among the most mTOR-kinase-inhibitor-sensitive sites in the compendium. It falls in 2 labs and 4 cell lines, while other LAMTOR1 sites are flat. This makes it a candidate direct mTORC1 site; S6K, SGK and other downstream kinases are not excluded."

## D2. RRAGC S15 (EETPLAG[S]YGAADSF)
**Reproduced and extended.**
- **CRC organoids (38047086, 2 h):** the mTORi gives −6.6 to −10.1 in **6/6 patient lines**. PI3Ki in the same organoids gives −0.7 (one unit), while EIF4B S422 falls to −5.2. This is the strongest evidence that the site is not on the AKT→AGC route.
- **BT474 (30901187):** dactolisib −5.9 and omipalisib −6.25.
- **RPE1 (32034161):** torkinib −2.2, but MEK inhibitors also give −3.0 (cobimetinib) and −2.3 (PD0325901), and AKT inhibitors −1.2 to −1.8. This study is nonspecific drift and should be discounted.
- The Cantley MTOR percentile is **100**, the best motif support of the four claims.
- Same-protein control RRAGC S95 is flat across all classes. The known mTORC1 feedback site T394 falls in single units (−8.7, −5.9).

**Weaknesses.**
- The rapalog/catalytic ratio of 0.48 rests on RPE1 alone: one study, a denominator of −2.2, and a noisy study. It carries little weight.
- No S6K-inhibitor, dose or CPTAC data exist.
- The organoid mTORi also drops LAMTOR1 S27 (−8.6) and T30 (−9.4). S27 is flat everywhere else, so the organoid condition has broad lysosomal or other effects.
- Effectively the evidence is 1.5 independent studies.

**Defensible wording.** "RRAGC S15 lies in the RagC N-terminal tail next to the established mTORC1 site S21, has an optimal mTOR motif, and is strongly mTOR-dependent in 6 patient organoid lines. The evidence is still provisional."

## D3. KPNA4 S60 (HEDICED[S]DIDGDYR): resolving the conflict with I26
**(i) The same experiments.** There are 4 CX-4945 conditions:
- Hijazi HL60, MCF7 and NTERA2: 10 µM, 1 h.
- HeLa (29214152): this is **Rusin et al. 2017, mitotically arrested HeLa**, 5 µM, about 45 min per the paper; the time is unannotated in the compendium.

| Condition | KPNA4 S60 z (rank) | PSP CK2 set median z (n) | In-vivo, LT≥3 CK2 sites, median | CK2-motif AUROC for z<−2 |
|---|---|---|---|---|
| HeLa CX-4945 | **−7.00 (63/12,455; 0.5th pct)** | **−1.33** (58); 41% below −2 vs 8% background | −0.24; 32% below −2 | **0.64 (best of 400 kinases)** |
| HL60 CX-4945 | −2.52 (416/8,111; 5.1th) | −0.01 (45) | +0.50 | 0.53 |
| NTERA2 CX-4945 | −4.12 (171/9,270; 1.8th) | +0.25 (48) | +0.50 | 0.49 |
| MCF7 CX-4945 | not measured | +0.41 (50) | +0.92 | 0.51 |
| TBCA × 3 lines | −0.04 / −1.64 / −0.16 | +0.04 to +0.21 | | 0.45–0.48 |

- **The I26 statement is correct for Hijazi at 1 h.** Known CK2 sites do not fall there as a set:
  - AKT1 S129 (the classic CX-4945 marker): HL60 +0.91, MCF7 −2.54, NTERA2 −0.29.
  - HSP90AB1 S226/S255, CTNNA1 S641, MCM2 S139 and SSB S366 are all near 0 or positive.
- **A fast subset does fall.**
  - Classic sites: NPM1 S125 −6.08 (NTERA2), EIF2S2 S2 −3.78 (NTERA2), OTUB1 S16 −5.46 (HL60), TOP2A S1525 −4.2 (HL60).
  - In HeLa, where CK2 inhibition is verified, the set falls as a whole.
- **CX-4945-specific module.** I took sites with z < −2 in both HL60 and NTERA2 that are down with ≤6% of the other ~58 Hijazi inhibitors. This gives 14 sites, including KPNA4 S60 (4%).
  - 64% have CK2 motif ≥90th percentile, vs 10% of background (**OR 16, p = 1.5e-6**).
  - CX-4945 off-target motifs are depleted: CLK1 OR 0.5, DYRK1A 0.2, HIPK2 0.4, CDK1 0.2.
  - Acidic module members: SF3A1 S451, ANLN S54, ZC3H8 S77, GPATCH2 S195, FXR1 S409, PAK2 T169, U2SURP S37, PHAX S368. FXR1 S409 is the paralog of I03's FXR2 S453.
  - The 6 module members with HeLa data have median z −6.0 there; 4/6 are below −2 vs 8% background.
  - The sites that CX-4945 lowers *non-specifically* are the S/T-P and RS-type ones (CLK/DYRK-like), and those are also lowered by other inhibitors.
- **Within Hijazi,** CX-4945 is the #1 of 61 inhibitors for KPNA4 S60 in NTERA2 and #3 of 55 in HL60.
- **Turnover.** KPNA4 S60 is a *quiet* site overall: |z| > 2 in 6.5% of 1,437 conditions (104 studies), vs a median of 12% for all sites and for PSP CK2 sites (19th–22nd percentile). Its CX-4945 response is therefore specific, not general volatility. It behaves like a member of the minority of CK2 sites with fast turnover; most CK2 sites are slow and do not move at 1 h, which is why I26's control fails.

**(ii) Off-targets** (Klaeger Kd_app for CX-4945):

| Target | Kd_app (nM) |
|---|---|
| CSNK2B / A1 / A2 | 0.15–1.0 |
| CLK3 | 10 |
| CLK1 | 14 |
| MAPKAPK3 | 26 |
| DYRK1A | 33 |
| STK10 | 96 |
| SLK | 220 |
| CDK5 | 514 |
| PLK1 | 672 |
| BMPR2 | 1,056 |
| PIM1 | 1,237 |
| CDK2 | 1,688 |

GSK3 and HIPK are absent below 3 µM in Klaeger (HIPK is listed only in Hijazi's own profile).
- **KPNA4 S60 motif percentiles for these kinases:** CLK1 6, CLK3 1, DYRK1A 4, HIPK2 1, PIM1 8, GSK3B 0, MAPKAPK3 5, CDK5 2, SLK 1, STK10 0. It does not fit any of them. It scores 99 for CSNK2A1/A2, and also 99 for CK1δ, GRKs and TGFβR/BMPR1, and 96 for PLK1. All are acidophilic.
- **Motif-set shifts** (top-2% motif sets) for CLK, DYRK, HIPK, PIM, GSK3 and MAPKAPK3 are ≈0 (|Δmedian z| ≤ 0.2) in every CX-4945 condition. Only CK2 shifts, and only in HeLa (−1.7).
- **Inhibitors of those kinases do not lower KPNA4 S60:**
  - CLK/DYRK/PIM/GSK3 inhibitors: BIM1, Go6983, abemaciclib and PF3758309 range from −1.4 to +3.4.
  - PLK1 inhibitors: BI2536 −0.48 / +0.68, BI4834 +0.42.
  - CHIR99021 and tideglusib: −1.5 to +0.3.
- **Unexplained outlier:** LY2090314 (GSK3i) in HL60 gives −6.17. It is not reproduced in NTERA2 or MCF7, by other GSK3 inhibitors, or by other module sites. I treat it as noise or an unknown off-target.

**(iii) Other CK2 inhibitors.** Only TBCA is present. It is weak and does not move any CK2 readout here, so it is uninformative rather than negative. There is no SGC-CK2-1, TBB, DMAT, quinalizarin, CX-5011, GO289, emodin or apigenin, and no genetic CSNK2 perturbation.

**(iv) Time and dose.** There are none: a single dose (10 µM; 5 µM in HeLa) and a single time point. decryptM has no CX-4945 curves. KPNA4 S60 is flat in all decryptM dose series (AZD8055, MK2206, MEKi, EGFRi).

**MAPK claim.** The refutation reproduces: 25 ERK-verified units show median +0.19 and none below −2, with ρ = −0.19 against the ERK activation loop. This addresses only acute (≤2 h) dependence; Hazawa used 24 h U0126.

**CPTAC.** Co-variation with the CK2 composite has meta-z 4.6. That is the 86th percentile of all sites but only the **58th among CK2-motif sites**, so it adds nothing beyond the acidic motif.

**Skeptic's view.** Menyhart et al. 2023 (inhibitor-resistant CSNK2A1 in U2OS) found that only about 15% of CX-4945-down sites at 4 h were truly CK2-dependent, vs >55% for SGC-CK2-1. A single chemotype, at 10 µM, without a rescue experiment, cannot prove CK2.

Counter-argument: the KPNA4 S60 module is acidic-motif-enriched 16-fold and depleted for every CX-4945 off-target motif, and its members fall hardest in the one experiment with verified global CK2 inhibition. Of the remaining acidophilic candidates, PLK1 inhibitors do nothing, and CK1, GRK and BMPR are not plausible CX-4945 targets at the observed potency. That leaves CK2 as the most parsimonious kinase.

**Verdict.** KPNA4 S60 is a **fast-turnover CK2-type site** with high confidence that the response is CX-4945-specific and medium confidence that it is CK2 itself rather than an off-target. No specific off-target is supported. It is not "unresolved" in the strong sense, but it awaits SGC-CK2-1 or an inhibitor-resistant CK2α rescue.

**Defensible wording.** "KPNA4 S60 belongs to a small set of acidic sites that fall within 1 h of CX-4945, specifically and in 3 cell lines. The set is 16-fold enriched for the CK2 motif and depleted of CX-4945 off-target motifs, and its members fall strongly in mitotic HeLa where global CK2 inhibition is verified. Acute MEK/ERK inhibition does not lower S60 (0/25 verified units). CK2 is the most likely kinase; this awaits a selective (SGC-CK2-1) or genetic confirmation."

## D4. ERK pathway → EPN1 S435 (AGEVPAR[S]PGAFDMS)
**Reproduced.**
- In 27 units with verified ERK-pathway inhibition (9 studies), the median z is −1.33 (9 units below −2); in units where the pathway is not down it is +0.24.
- Against 2,164 S/T-P sites this is the **1.2th percentile**. Only 2.2% of S/T-P sites have a median ≤ −1.
- Spearman ρ with the ERK pT/pY composite is 0.57 (n = 52), the **99.4th percentile** of the S/T-P null (99th percentile = 0.51). LOSO gives ρ 0.36–0.62.
- **The median depends on decryptM.** Dropping decryptM (36926954) moves the median-in-down-units to −0.53; dropping Hijazi moves it to −1.96.
- Single strongest units: trametinib in WM239A −10.2, selumetinib in HL60 −8.7, ulixertinib in HL60 −7.9, PD184352 in A375 −5.5.
- **Inconsistent units:** no drop under MCF7 trametinib/selumetinib despite full ERK loss, PC9 EGFR inhibitors, A549 selumetinib, or ravoxertinib.
- **CDK inhibitors** have no effect (+0.12, 6 units).
- **Same-protein controls:** T470, S447 and T460 are flat. S454 (also S/T-P) co-responds (ρ 0.64). That is consistent with ERK acting on two S-P sites in the same disordered C-terminus, but it means S454 is not a negative control.

**Corrections to V-I03.**
- The "potency matches ERK" claim relies on refametinib pEC50 7.80. That fit is **not significant** (Δ −0.59), and the PD0325901 A549 fit is flat (Δ −0.2).
- Only the lapatinib fits are significant: MDA-MB-175-VII 7.05 and BT474 5.85.
- The dose evidence is therefore weak.

**CPTAC.** Co-variation with ERK has meta-z 6.5 (96th percentile among S-P sites), but with the CDK composite it is 8.5 (98th percentile).

**Defensible wording.** "EPN1 S435 is a reproducible, acute MEK/ERK-pathway-dependent site: 99th percentile of S/T-P sites, 9 studies. It is a candidate direct ERK1/2 substrate, but direct phosphorylation and the exact proline-directed kinase are untested."

## Literature check (2026-10-29, brief)
- **KPNA4 S60:** the only functional paper is Hazawa et al., Oncogene 2019 ("MAPK-dependent", from 24 h U0126 and a pan-pSer pulldown). I found no CK2 report.
- **LAMTOR1 S97/S98:** no kinase or functional report.
- **RRAGC S15:** no report. Yang et al., EMBO J 2018 (30552228) established mTORC1 → RagC S21 autoregulation, which supports plausibility.
- **EPN1 S435:** no report. Only mitotic and Cdk1 epsin phosphorylation is known.
- **Relevant caution:** Menyhart et al., Curr Res Chem Biol 2023 found that CX-4945-down sites are mostly CK2-independent, which motivates SGC-CK2-1 validation.
- I could not open full texts (PMC was blocked by a captcha), so I could not confirm whether KPNA4 S60 appears in CK2-rescue datasets.

## Recommended next experiments
1. KPNA4 S60: SGC-CK2-1 time course (5–60 min) with the inhibitor-resistant CSNK2A1 (Menyhart line), with SF3A1 S451 and ANLN S54 as module controls.
2. LAMTOR1: targeted PRM to resolve S96/S97/S98; rapamycin vs Torin1; S6K1/2 double knockout or PF-4708671 with verified RPS6 S240/244 loss.
3. RRAGC S15: in vitro mTORC1 kinase assay on an RagC N-tail peptide (S15A/S21A).
4. EPN1 S435: in vitro ERK2 assay; ERK inhibitor panel with RSK S363 readout.

## Files
- `tables/kpna4_ck2i_rank.tsv`, `ck2_psp_sites_in_ck2i.tsv`, `ck2_set_summary.tsv`, `motif_set_shift_ck2i.tsv`, `motif_auroc_decreased_ck2i.tsv`
- `tables/kpna4_s60_all_conditions.tsv`, `kpna4_offtarget_inhibitors.tsv`, `cx4945_consistent_down_hijazi.tsv`, `cx4945_specific_module.tsv`
- `tables/class_effects_loso.tsv`, `units_mtor_candidates.tsv`, `s6ki_efficacy.tsv`, `decryptM_panel.tsv`, `epn1_erk_units.tsv`, `cptac_null.tsv`
