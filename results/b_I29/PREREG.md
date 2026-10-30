# I29 PREREG — ERK-phosphorylation "ceiling" in EGFR-mutant LUAD
Written 2026-10-29 BEFORE any new outcome was computed. Already known (PI preliminary, b_I27): EGFR-mut vs WTall stage_ERK −0.78 (padj), ERKSUB −0.68, ETV4 +0.73, ETV5 +0.45, SPRED1 +0.64, DUSP6 +0.50, DUSP4 −0.58. Inspected before writing: only file layouts, column names, group counts, the existence of Broad Spectrum-Mill multi-site phosphopeptide IDs for MAPK1/MAPK3 (T185t, Y187y, T185tY187y; T202t, Y204y, T202tY204y) and the compendium condition list (PC9 gefitinib dose series; H1975/H3255 in PMID 25404012). No values were looked at.

## Groups (CPTAC LUAD, n=110 tumours)
- EGFR-mut: any non-synonymous EGFR mutation (n=39; PI definition). Clean EGFR: ex19del or L858R (incl. multi containing them), no KRAS/NRAS/BRAF/NF1 co-mutation.
- RAS-path: KRAS/NRAS/BRAF/NF1 mutant without EGFR mutation (n=42).
- WTall: none of the above (n=29; reference). Clean WTall: also excluding ALK/RET/ROS1 fusions, EGFR-AMP and ERBB2 in-frame (cBioPortal PUTATIVE_DRIVER_MUTATION / *_FUSION), and EGFR GISTIC = 2.

## Outcomes
- pERK: stage_ERK (BCM site-level, tc_padj) [primary phospho]; Broad Spectrum-Mill peptide-level MAPK1/MAPK3 mono-pT, mono-pY and doubly pTpY separately, protein-adjusted by residualising on Broad ERK protein (same gene) [key secondary]; pTpY fraction index = pTpY − mean(mono).
- ERKSUB, RSK (stage_RSK), b_I10 footprints fp_KI_MEK (sign-flipped = MAPK output), fp_EGF, fp_MAPKaxis.
- ETV composite (RNA) = mean of within-cohort z of ETV4, ETV5, SPRED1, DUSP6 [primary transcriptional]; individual genes plus ETV1, SPRY2, SPRY4, DUSP4, PHLDA1, EGR1, FOS (IEG; ischaemia-sensitive), and PROGENy MAPK.
- Proteins (BCM/Broad) where measured: MAPK1, MAPK3, DUSP4, DUSP6, SPRY2/4, SPRED1/2, ETV4/5.

## Model
OLS outcome ~ group + purity + E2F (prolif) [base], HC3 95% CI, Freedman–Lane permutation p (5,000) for the two primaries; HC3 p otherwise. Effects in within-cohort SD.

## Gate G — cold-ischaemia / collection-site confound (coordinator-mandated; run first)
- ISCH = audit_I10 ischaemia footprint (top-200 Mertins 2014 cold-ischaemia sites, `scripts/audit_I10/common.py::ischaemia_set`), scored as signed mean of standardised input_tc_padj LUAD tumour sites (`score`). ISCH_noERK = same after removing sites on MAPK1/3, RPS6KA1/3, PSP ERK1/2 substrates and any site in b_I10 KI_MEK / KI_ERK / KI_RAF / EGF / GF_all / MAPKaxis footprints (to avoid adjusting away the outcome itself).
- G1: ISCH, ISCH_noERK by group (base model; Kruskal–Wallis as well).
- G2: re-run pERK, ERKSUB, ETV composite, ETV4, ETV5 with + ISCH_noERK; + ISCH; + collection-centre proxy (COUNTRY_OF_ORIGIN: Vietnam / USA / China / other) ; + both; within-Vietnam only; + TMT plex fixed effects (sensitivity). Ischaemia times / tissue source site are not in local tables; a GDC/PDC biospecimen query (time_between_clamping_and_freezing, tissue source site code) will be attempted and used if obtained.
- Verdict: dissociation SURVIVES if in the fully adjusted model (base + ISCH_noERK + country) pERK β<0 with p<.05 AND ETV composite β>0 with p<.05; ATTENUATED if same signs but either p≥.05; FAILS if a sign flips or |β_pERK| falls >70%.

## Hypotheses (one-sided where direction stated)
- H1 (primary): EGFR-mut < WTall on pERK; H1b: on doubly-phosphorylated pTpY of MAPK1 and MAPK3 (Broad, protein-adjusted). Mono-pT / mono-pY: two-sided (DUSP6 removes both; a PTP-type cap would leave mono-pT up).
- H2 (primary): EGFR-mut > WTall on ETV composite.
- H3: RAS-path vs WTall pERK not lower (β ≥ −0.2) — the "KRAS behaves differently" claim; also EGFR vs RAS-path direct contrast.
- H4 robustness: H1 β<0, p<.05 survives each of: + EGFR allele (ex19 vs L858R within EGFR), + smoking, + stage, + sex, + TP53, clean EGFR vs clean WTall, tc instead of tc_padj.
- H5 (within EGFR-mut, n≈39): DUSP6 RNA (and protein if measured) correlates negatively with pERK (Spearman, partial on purity) — informative; and positively with ETV4/5 (expected trivially since DUSP6 is an ERK target; reported, not counted as support).
- H6 LSCC & other CPTAC cohorts with EGFR mut/amp: descriptive, same direction expected only for mutation (amp-driven tumours e.g. GBM are not predicted).
- H7 replication (if obtainable within 30 min): EGFR-mut vs EGFR-WT (or driver-negative) lower pERK, one-sided p<.05, and higher ETV4/5 RNA if RNA available.
- H8 compendium: in EGFR-mutant lines under EGFR-TKI, MAPK1/3 activation-loop sites fall with dose (Spearman ρ<0) — "low-but-necessary ERK". Basal level comparisons across lines are not possible from log2FC data (declared untestable).

## Pre-declared limits
n per group 29–42 (MDE ≈ 0.7 SD). Country of origin is collinear with EGFR status (Asian enrichment) — reported as such. Peptide-level data are from a different pipeline (Broad) than BCM site-level and may have more missingness.

## Amendment A1 (before any outcome computed, 2026-10-29)
- GDC biospecimen query succeeded: `time_between_excision_and_freezing` (cold ischaemia, min; median over tumour aliquots) for 106/110 LUAD cases, `time_between_clamping_and_freezing` (warm+cold) for 70. Added to Gate G: (i) excision→freezing time by group; (ii) correlation of time with pERK and with the ISCH score (validates the footprint as an ischaemia proxy); (iii) primary contrasts additionally adjusted for log(excision time). Tissue-source-site codes are not exposed by GDC for CPTAC-3; COUNTRY_OF_ORIGIN remains the centre proxy.
- Coordinator supplied val_E ischaemia scores (`results/val_E/ischaemia_scores.tsv`: isch, isch_noPath = pathway-protein sites removed). Used as given alongside my stricter ISCH_noERK (also removes ERK substrates/MAPK-footprint sites). Fully adjusted model = base + ISCH_noERK + country + log(excision time).

## Amendment A2 — replication cohort (declared after Gate G / CPTAC robustness were run, BEFORE any replication data value was inspected)
- Cohort: Chen et al. 2020 Cell (Taiwan, TW LUAD, n≈103 patients), Table S1 (Cell mmc1.xlsx; S1A clinical incl. EGFR_Status; S1C SNV; S1D RNA log2 T/N; S1E proteome log2 T/N; S1F phosphopeptide; S1G phosphosite log2 T/N). Secondary: PDC000220 CDAP phosphosite.tmt10 (tumour vs pooled reference) if tumour channels can be identified.
- Groups: EGFR-mut = EGFR_Status not WT; reference = EGFR_Status WT and no KRAS/NRAS/BRAF/HRAS/NF1 non-synonymous SNV ("WTall-like"; fusions unknown); RAS-path = KRAS/NRAS/BRAF/NF1 mutant, EGFR WT. Primary set: adenocarcinoma histology only (S1A Histology Type == ADC); all tumours as sensitivity.
- Outcomes: MAPK3 T202/Y204 & MAPK1 T185/Y187 (site-level S1G; pTpY doubly-phosphorylated peptide from S1F if present), pERK = mean z of available activation-loop sites; protein-adjusted version = residual on S1E MAPK1/3; ETV composite (ETV4, ETV5, SPRED1, DUSP6 RNA T/N); DUSP6, DUSP4 RNA/protein.
- Model: OLS outcome ~ group (+ stage, + smoking as sensitivity); HC3. Replication = pERK β<0 one-sided p<.05 AND ETV composite β>0 one-sided p<.05. Caveat declared: T/N ratios differ from CPTAC's tumour-vs-reference scale; no purity covariate.
