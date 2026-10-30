# I27 PREREG - allele-specific signalling outputs (CPTAC)
Written 2026-10-29 BEFORE computing any outcome. Only allele counts per cohort and feature names were inspected beforehand.

## Data / definitions
- Alleles from BCM MAF (mutations.parquet). Transcripts are canonical: KRAS ENST00000256078 (4B, P01116), PIK3CA ENST00000263967 (P42336), BRAF ENST00000646891 (MANE, P15056, V600 present), EGFR ENST00000275493 (P00533, L858R/E746_A750del present). Codon numbering therefore matches UniProt.
- Allele groups: KRAS G12D / G12V / G12C / G12R / G13D / Q61H. PIK3CA helical = E542K, E545K, E545A, Q546x; kinase = H1047R/L/Y; double = >=2 distinct PIK3CA missense (excluded from single groups). BRAF class1 = V600E/K; class2 = G469A/R/V, K601E, L597x, in-frame b3-aC indels; class3 = D594x, G466x, N581x, G596x. EGFR ex19del = in-frame del/delins at codons 745-753; L858R separate.
- WT = no non-synonymous mutation in that gene. Tumours carrying two different alleles of the same gene are dropped from allele-vs-allele contrasts.
- Readouts, primary variant tc_padj (protein-adjusted), tc as sensitivity:
  (a) I24 stage composites RTK, EGFRpY, ADAPT, ERK, RSK, ERKSUB, AKTSUB, MTORC1, S6;
  (b) kinase activity corev3_z (family-level);
  (c) I10 footprint scores (nz, 49 groups);
  (d) pre-specified sites: AKT1S1/PRAS40 T246 (AKT-specific) and NDRG1 T328+T366 mean (SGK-specific).
- Each outcome z-scored within cohort across all tumours.

## Model
Per cohort OLS: outcome ~ group + purity + E2F (+ subtype where stated). Reported: effect in within-cohort SD, HC3 95% CI, Freedman-Lane permutation p (5000 perms), MDE = 2.8*SE (80% power, two-sided a=.05). Across cohorts: inverse-variance fixed-effect and DerSimonian-Laird random-effects meta.

Verdict rules (pre-declared):
- Confirmed: predicted direction and p<.05 (meta, or primary cohort).
- Contradicted: opposite direction and p<.05.
- Not confirmed: p>=.05 and MDE <= 0.8 SD.
- Underpowered: p>=.05 and MDE > 0.8 SD.

## Predictions (sources verified by web search; corrections noted)
P1. Hobbs et al. 2020 Cancer Discovery 10:104 (NOT Nature, as in the brief): KRAS G12R interacts poorly with p110a.
 - Primary: PDAC, G12R vs G12D+G12V, outcome AKTSUB, additionally adjusted for ERK stage ("matched MAPK output"). Predicted negative.
 - Secondary: MTORC1, S6, AKT-family activity, PRAS40 T246 lower; ERK/ERKSUB not different.
P2. Ihle et al. 2012 JNCI 104:228 (verified): G12C and G12V lines show Ral activation and reduced GF-dependent AKT; G12D shows PI3K and MEK activation. Source groups G12V WITH G12C, so registered as:
 - Primary: G12C < G12D on AKTSUB; LUAD primary, meta with COAD/UCEC where n>=3.
 - Secondary: G12C vs G12V null; G12V < G12D on AKTSUB (meta PDAC, COAD, UCEC, LUAD). ERK not different.
 - RAL arm untestable (no RalGEF/Ral phospho readout) - declared untestable in advance.
P3. Zhao & Vogt 2008 PNAS 105:2652 (verified): helical mutants require RAS-GTP, kinase-domain mutant requires p85. AKT-vs-SGK axis from Vasudevan 2009 Cancer Cell 16:21 (low AKT, PDK1/SGK3 dependence in PIK3CA mutants; NOT allele-resolved). Vasan 2019 Science concerns double cis mutants, not helical vs kinase. No source states a direction for helical vs kinase, so P3 is registered TWO-SIDED:
 - Primary: helical vs kinase on (SGK - AKT) contrast = NDRG1 composite minus PRAS40 T246, and on AKTSUB and MTORC1; adjusted for PAM50 (BRCA) and MSI/genomic subtype (UCEC, COAD). Meta BRCA, UCEC, COAD, HNSCC.
 - P3b (Vasan 2019): double PIK3CA > single on AKTSUB (one-sided, pooled).
P4. Yao et al. 2015 Cancer Cell / 2017 Nature 548:234 (verified): class 1 (V600) suppresses RAS-GTP via ERK feedback; class 3 is RAS/RTK-dependent.
 - Primary: COAD, BRAF V600E vs KRAS/NRAS-mutant tumours (the RAS-dependent comparator, since class 2/3 n ~ 1-3 per cohort): V600E predicted LOWER RTK pY and ADAPT pY.
 - Secondary: V600E vs triple-WT in COAD, predicted lower; class 2/3 pooled pan-cancer vs V600E, predicted higher RTK/ADAPT. Declared in advance as likely underpowered.
P5. EGFR ex19del vs L858R, LUAD (n 15 vs 15). Literature search found NO primary source reporting a specific adaptor/pathway direction (Chen 2006 Oncogene tested other mutants; clinical PFS differences only). P5 is therefore registered as DIRECTIONLESS/two-sided on EGFRpY, RTK, ADAPT, ERK, AKTSUB, MTORC1, and EGFR/ERBB family activity, and the brief's framing is noted as not literature-anchored.

## Exploratory
All outcomes (stages + kinase activities + footprints) x all allele pairs with n>=8 per arm in the same cohort, plus each allele vs gene-WT. BH FDR within the exploratory family. Reported with n, effect, CI, q.

## Pre-declared limits
Kinase activity is family-level; receptor pY tracks receptor abundance (so tc_padj is primary). Comparisons are small-n; most pairwise tests will have MDE around 0.8-1.4 SD and verdicts will reflect that.
