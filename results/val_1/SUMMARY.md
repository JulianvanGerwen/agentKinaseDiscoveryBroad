# val_1: cross-validation of V1-V3 (compendium / CPTAC), 2026-10-29

Scripts: /workspace/scripts/val_1/ (v1.py, v2.py, v2b.py, rescore_masked.py, figs.py). Tables: tables/. Figures: figures/V1_reporters_compendium.png, V2_hippo_meta.png, V3_pak4_s148.png.

## V1. Tumour-derived reporter sites in perturbation data: MOSTLY UNTESTABLE; weakly supported where testable
Compendium has almost no on-target inhibitor conditions for these sites (tables/V1_inhibitor_effect.tsv). Effect = median log2FC of group minus all other qc-ok conditions where the site is measured; AUROC = P(group lower than other) (stimulus: higher); within-study = same comparison inside studies containing both.

| site | group (n cond / n studies) | effect | AUROC | p |
|---|---|---|---|---|
| ABL1 Y185 | imatinib/dasatinib/crizotinib (26/3; 22 are decryptM) vs 18 others | -0.04 | 0.58 | 0.40 |
| CDK13 T588 | 0 CDK-inhibitor conditions among 514 (THZ531/dinaciclib not measured at this site) | NA | NA | untestable |
| PTK2 S840 | crizotinib/ceritinib (3/2; ALK drugs, off-target FAK, no FAK inhibitor) | -0.75 | 0.88 | 0.018 (n=3) |
| PTK2 Y576 | same drugs (7/4) | +0.06 | 0.38 | 0.27 (wrong direction) |
| PAK2 S197 | BIM1 only (3/2; not a PAK inhibitor), no PF-3758309/FRAX measured | -0.80 | 0.68 | 0.29 (untestable) |
| IKBKB S697 | SC514/TPCA-1/MRT (5/2) vs 54 | -2.5 | 0.81 | 0.021 |
| IKBKB S697 | TNF (29/3) vs 54 | +1.39 | 0.75 | <0.001 |
| IKBKB S675 | IKKi (2/1): n too small; TNF (5/2) went DOWN (-0.84, AUROC 0.15) | | | 0.007 (wrong direction) |
| MAPK1 Y187 (control) | ERKi (19/4) / MEKi (58/7) | -0.69 / -1.43 | 0.72 / 0.77 | 0.001 / <1e-4 |
| MAPK3 Y204 (control) | ERKi (16/3) / MEKi (53/7) | -0.50 / -1.06 | 0.70 / 0.77 | 0.005 / <1e-4 |

Activity coupling across all conditions (tables/V1_activity_corr.tsv). Important: the compendium core-v3 and PSP scores contain each reporter site as a substrate of its own kinase (autophosphorylation edges; checked in the edge files), so the stored scores are circular. I re-scored core-v3 with all sites of the 7 kinase proteins masked ("core_masked") and also computed a PSP mean-z excluding own-protein sites ("loso", needs >=5 substrates).

| site | n cond (studies) | rho core_masked | rho loso-PSP (n) |
|---|---|---|---|
| ABL1 Y185 | 44 (5) | 0.15, p=0.33 | 0.25 (31), p=0.18 |
| CDK13 T588 | 512 (20) | 0.42, p<1e-3 | -0.13 (103, 3 studies), p=0.20 |
| PTK2 S840 | 384 (29) | 0.25, p<1e-3 | not scorable (<5 non-self substrates) |
| PTK2 Y576 | 449 (33) | 0.22, p<1e-3 | not scorable |
| PAK2 S197 | 576 (66) | 0.22, p<1e-3 | 0.33 (65), p=0.008 |
| IKBKB S675 | 212 (17) | -0.11, p=0.10 | -0.09 (103), p=0.35 |
| IKBKB S697 | 86 (12) | 0.25, p=0.02 | 0.52 (24), p=0.009 |
| MAPK1 Y187 (control) | 1484 (93) | 0.47 | 0.28 |
| MAPK3 Y204 (control) | 1337 (85) | 0.32 | 0.29 |

Within-study-centred values are similar (tables). Caveats: CDK13 core-v3 score leans on CDK-family borrowing, so the 0.42 may reflect general CDK/proliferation activity, and the substrate-based PSP check (3 studies) does not support it. All kinases have far fewer conditions than the MAPK controls.
**Verdict:** IKBKB S697 supported (down under IKKi n=5, up under TNF, coupling 0.25-0.52); PAK2 S197 and PTK2 S840/Y576 show modest positive activity coupling (rho 0.2-0.3, about half the MAPK control strength) but inhibitor tests are untestable (off-target or no relevant drugs). ABL1 Y185 not supported (no effect under imatinib/dasatinib, no coupling). IKBKB S675 not supported (wrong direction under TNF). CDK13 T588: inhibitor test untestable, coupling ambiguous. Overall: partial, mostly untestable; only MAPK1/3 and IKBKB S697 are validated.

## V2. MET/ALK inhibition to Hippo/NDR activation: NOT SUPPORTED (opposite sign in tumours; YAP1 S127 untestable)
10 CPTAC cohorts (tc input), partial Spearman adjusted for purity (WES, ESTIMATE fallback) and proliferation (Hallmark G2M+E2F), random-effects meta (tables/V2_meta.tsv). Hippo = mean z of LATS1/LATS2/STK38/STK38L core-v3. "noYAP1" = re-scored with all YAP1 sites masked, and MET/ALK/MAPK1/3/Hippo-kinase own-protein sites masked. Hypothesis predicts r<0.

| predictor -> Hippo score | meta r [95% CI] | cohorts r<0 |
|---|---|---|
| MET activity -> full | -0.05 [-0.15, 0.05], p=0.33 | 6/10 |
| MET activity -> noYAP1 | -0.08 [-0.16, 0.01], p=0.08 (STK38 alone -0.09, p=0.04, nominal) | 7/10 |
| ALK activity -> full / noYAP1 | +0.23 / +0.14 (p<1e-3 / 0.03) | 2/10, 3/10 |
| RTK/MAPK axis (EGFR, MET, ALK, ERBB2, MAPK1/3) -> full / noYAP1 | +0.23 / +0.19 (p<1e-3 / 1e-3) | 2/10, 3/10 |
| MET protein -> full / noYAP1 | +0.02 / +0.06 (ns) | |
| MET-high (top decile protein) vs rest, Cohen d | +0.14 [-0.07, 0.35] / +0.20 [-0.04, 0.44] (ns) | 2/10, 3/10 |

Removing YAP1 sites barely changes the picture (the Hippo score stays at the same sign and about the same size). YAP1 S127 (protein-adjusted) is essentially unmeasured: tumour-centred phospho present in only GBM (n=34) and OV (n=35); protein-adjusted only in GBM. In GBM, MET activity vs YAP1 S127 r=+0.38 (p=0.035), wrong direction, single cohort, no meta possible.
Caveats: MET/ALK scores rest on a median of 9/8 substrates (largely predicted); the tumour test is cross-sectional, while the claim is about acute inhibition, so this is an indirect test. Positive RTK/MAPK-Hippo coupling may reflect shared substrates/proliferation/stromal structure. A negative association is not seen, so the claim gets no tumour support.

## V3. PAK4 S148 up under PF-3758309: UNTESTABLE for replication (single study), plausibly nonspecific
- S148 is measured in 734 conditions / 55 studies, but PF-3758309 exists only in Hijazi (PMID 31959955; HL60 +2.45, MCF7 +3.06, NTERA2 +2.59; within-cell-line percentile 97-98%, robust z 1.7-2.9). No KPT-9274 or LCH-7749944 anywhere; the only other PAK-directed drug is FRAX486 (group I PAK, same study), which moves S148 DOWN (-0.6 to -1.1). Other PAK inhibitors not present.
- Not specific: in Hijazi alone 4/165 conditions reach >=2.45 (ulixertinib MCF7 +3.5, abemaciclib NTERA2 +2.2) and the site is high in mitosis vs S/control (+2.5, +2.8; 2 studies) and under vemurafenib(+gefitinib) 24 h (+2.2 to +2.6). Non-Hijazi inhibitor conditions: 361, median ~0, SD 0.39. This suggests a cell-cycle/MAPK-linked response (kinase-substrate prediction: CDK1, p=0.11) that PF-3758309 may produce via arrest or off-targets (it hits many kinases) rather than a PAK4 self-site effect. Not tested further.
- Other PAK4 sites under PF-3758309 (Hijazi): S474 (activation loop, PSP regulatory) -1.4/-0.6/-1.0, consistent with on-target inhibition; S181 ~0; S104 inconsistent.
- PSP: S148 is in Phosphorylation_site_dataset (18 high-throughput records, no low-throughput, no kinase) but absent from Regulatory_sites (no known function). Literature search (PubMed-style web query on Ser148/PF-3758309) found nothing on this site; the PNAS 2010 paper (Murray, PMID 20439741) reads out GEF-H1 phosphorylation; PF-3758309 is flagged as promiscuous (13/146 kinases hit).
**Verdict:** untestable by replication, and the single-study signal looks nonspecific; keep as a hypothesis only if a PAK4-selective inhibitor (KPT-9274) or PAK4 KO experiment is run.
