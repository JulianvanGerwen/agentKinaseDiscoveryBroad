# RESEARCH TREE — Discovery in kinase signaling (discoveryBroad2)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-29 FINAL (T+~18h) — claims frozen, REPORT.md written
Claims were frozen after the last audits (A-F, A-H, V-B) and REPORT.md was written from audited wording only. There are four key findings. **K1:** KPNA4 S60 is a CK2α-dependent, fast-turnover site; this corrects a published MAPK attribution and is confirmed by independent public chemical-genetic rescue data. **K2:** signalling programmes learned from cell-line perturbations are preserved in human normal and tumour tissue (10 tumour and 6 NAT cohorts; 91% sign agreement), not explained by measured ischaemia, and tumours specifically gain DNA-damage/replication-stress coherence. **K3:** primary-direction phosphoproteome responses are portable across labs; reverse-direction responses are study-specific (a methods caution). **K4:** in tumours, receptor-abundance-driven pY is buffered from ERK while acute perturbations propagate, and lesions shift cascade offsets, not gains. Caveated observations: candidate mTORC1 sites on Ragulator/Rag (LAMTOR1 S97/98, RRAGC S15), the turnover atlas, the ATP/Cheng–Prusoff consistency check, the PSP audit, and activity tracking activating sites in tumours.

About 24 branches and 7 audits ran; most were negative and are documented. Audits changed the wording of every key finding, and two candidates (G: ERK ceiling; immune–PTP mimicry) were dropped after replication or audit. An over-adjusting ischaemia proxy was caught and withdrawn using measured GDC times. An independent fact-check of REPORT.md and the final retro are running (docs/retro/FINAL_factcheck_retro.md).

### 2026-10-29 Milestone M2 (T+~14h)
Since M1, 14 more branches, validations and audits have run. Candidate findings with independent audits now span four data types. **A** (audited 3×, replicated in I10 as classification and in tumours): primary-direction phosphoproteome responses reproduce across labs and cell lines; reverse-direction responses reproduce only within a study's culture. No mechanism found; cell-cycle, stress, Hippo, normalisation, baseline expression and drug resistance were all ruled out. **B** (audited, qualified): cell-line perturbation footprints are coherent units in all 10 CPTAC cohorts at ~4× PSP sets; the MAPK axis is half shared with a cold-ischaemia signature, while the DNA-damage footprint is clean. **D** (validated + audited): named edges. KPNA4 S60 is a fast-turnover CK2 site, correcting a published MAPK attribution, inside a 14-site fast CK2 set; LAMTOR1 S97/98 and RRAGC S15 are candidate direct mTORC1 sites on the lysosomal docking machinery; EPN1 S435 is ERK-pathway. **E** (audited, reworded, ischaemia-robust): in tumours, abundance-driven receptor pY does not reach ERK (more receptor at fixed pY → lower ERK, 6/7 cohorts), while acute receptor perturbations in cell lines propagate. **F**: a perturbation audit of PSP finds 27 corrections, mostly AKT/CDK/PKC-annotated sites that behave as ERK/RSK sites, several literature-confirmed. Secondary: an ATP-competition consistency check (I17); lesions shift cascade offsets, not gains (powered negative, I24).

New lead **G** (I29 running): EGFR-mutant LUAD has lower ERK phosphorylation but higher MAPK transcriptional output (ETV4/5, SPRED1) than driver-negative LUAD, with low DUSP4. This fits the ERK-toxicity ceiling model and is ischaemia-robust; replication in an independent LUAD cohort is being attempted. Negatives (dark kinases, phosphatase natural experiments, CNV-IV causal edges, allele predictions, reverse layer vs resistance or mechanism) are documented. Remaining work: I29, Q8 (tumour vs normal wiring), I28 (turnover atlas); then the final report with a ranked findings table, a final retro and lesson pruning.

### 2026-10-29 Milestone M1 (T+~6h)
Infrastructure (literature, compendium, CPTAC) and seven discovery/validation branches are complete. The leading result is **Finding A (candidate key, audited)**: phosphoproteome responses to perturbations have two layers. The primary-direction layer (decreases under inhibitors, increases under stimuli) replicates across cell lines and studies, more so with effect size, is enriched for motif-matched direct substrates, and is ~2× more potent in dose series. The reverse-direction layer replicates within a cell line but not across cell lines; its kinase-motif signature (NEK, IRAK4, STE20/MAP4K, MAP3K, MARK/SIK) is cell-line dependent. The stimulation mirror test is decisive; the inhibitor arm alone leans on MEK, and timing does NOT explain the layers (I05 negative). Secondary results: downstream phosphosite potency shifts imply ~2 mM intracellular ATP via Cheng–Prusoff (I17, n=5 drugs); kinase activity in tumours follows activating-site phosphorylation, not abundance (I08, partly known).

Negatives so far: cis-CNV IV causal edges (I20; kinase abundance not limiting), cross-drug off-target assignment (scooped by decryptM 2.0), MET/ALK→Hippo (not in tumours), PAK4 S148 and most tumour reporter sites (failed cross-validation). Running: I10/I06 phospho-CMap + modules, I01 dark kinases, I22 cell-intrinsic reverse repertoire (explains Finding A's context dependence?), I23 reverse layer vs drug resistance (translational relevance of Finding A). Next decisions: if I22/I23 positive, consolidate Finding A with an independent audit and figure; seek at least two further independent discoveries (tumour-side from I10, dark kinases from I01); the decryptM 2.0 access request is pending with the user.

### 2026-10-29 T+0 — kickoff
Project started from PROMPT.md. Located the two method repos at
`<HACKATHON_ROOT>/{kinaseSubMain,kinaseActivityMain}` (read-only):
kinaseSubMain = kinase–substrate predictor (released `results/release/site_top5.tsv.gz`, `kinase_top200.tsv.gz`; 503 kinases);
kinaseActivityMain = core-v3 activity inference (`mvp/`, Bayesian + family prior + signatures, 503 kinases, ~7 ms/condition).
Local data: compendium of 1,937 conditions (data/), PSP, Klaeger, Hijazi inhibitor targets; plus in kinaseSubMain/external:
CPTAC pan-cancer (10 cohorts, v1.2 + PDC harmonized), DepMap, GTEx, HPA, interactomes (BioPlex, HuRI, OpenCell, BioGRID, STRING), Cantley atlas, Ochoa 2020 functional scores.

Wave 1 launched: literature review, compendium infrastructure, CPTAC infrastructure. Next: idea generation + ranking from the literature and data profile, then wave 2 of discovery branches.

## Decision log
- D0 (T+0): Resource cap = a QUARTER node (≤16 cores, ≤94 GB RAM, ≤1 GPU), stricter than PROMPT's "half": the user set this on 2026-10-28 (memory) because other agents share the node. Split across concurrent agents.
- D1 (T+0): Method repos are used as tools (read-only); outputs go to /workspace. Known tool limits carried forward: Tyr-kinase substrate predictions unreliable; activity inference cannot separate kinases within a family → claims at kinase level need family-level sanity checks.

- D2 (T+~1h): Wave 2 = I17/18/21 dose pharmacology (decryptM×Kinobeads), I19 feedback & self-hyperphosphorylation atlas, I10/I06 phospho-CMap + modules→CPTAC, I20 MR-style causal edges in CPTAC. Dropped I14 (circular), downgraded CPTAC driver/survival claims (largely published).

- D3 (T+~1.2h): Library mismatch (compendium = evaluated curated_union; CPTAC = MVP PSP+OmniPath recipe, 91% curated overlap) accepted: cross-dataset work will operate at site level or re-score with one library when needed.

- D4 (T+~3h): I19.1 promoted to key-finding candidate after audit (mirror test). Follow-up B-I05: test whether primary/portable responses are also early and motif-direct (time series), unifying I19 with the temporal-ordering gap (lit Q15).

- D5 (T+~5h): Consolidating I19+I05+I17 into candidate Finding A, a two-layer architecture of phosphoproteome responses: a portable, primary-direction, motif-direct, potent layer and a cell-line-specific, reverse-direction, less potent layer. Timing does NOT separate them. Next: explain the cell-line dependence of the reverse layer from baseline cell-line profiles (B-I22).

- D6 (M1 retro, docs/retro/M1_retro.md): (1) Auditors, not producing agents, assign C scores; [!] only with audited wording. (2) Rebalance away from Finding A toward named biological discoveries: launch B-I03 (new kinase–substrate edges by integrating co-regulation, motif, dose potency, interactome; held-out known edges for precision) and B-I13 (phosphatases) as cores free. (3) Gate I10: tumour-side claims must survive subtype/proliferation/purity adjustment, else drop. (4) ~2 mM ATP is a consistency check, not a measurement.

- D7 (T+~7h): Finding A reframed after I22: 'primary-direction responses reproduce across labs/cell lines; reverse-direction responses reproduce only within a study's culture' — a caution for interpreting feedback/compensation from single studies. Not 'cell-intrinsic compensation'. Finding A stays a candidate; I23 (resistance) decides whether the reverse layer carries biology.

- D8 (M2 retro, docs/retro/M2_retro.md): Finding A CLOSED (methods/reproducibility caution). Before write-up: independent audit of F + public-data CK2-selective check of KPNA4 S60 (A-F); G needs replication + independent audit. Then freeze claims (lesson 13).

- D9 (FINAL, corrected after fact-check docs/retro/FINAL_factcheck_retro.md): Tier 1: K1 KPNA4 S60 CK2 (audit_D 3.5/3/3.5, raised to high confidence by audit_F public rescue). Tier 2: K2 preserved perturbation wiring + tumour DDR gain (audit_I10 qualified; Q8 unaudited), K3 primary vs reverse replication asymmetry (audit_I19 3/4/3, methodological), K4 abundance buffering + offset-not-gain (audit_I24 qualified; PI score 3/3/4). G dropped (pERK half failed CPTAC Confirmatory; Taiwan underpowered).

## Ranking scheme (applied to every result)
Each candidate finding scored 1–5 on: **Compelling** (effect size, replication across independent datasets/cohorts, robustness to method choice, beats null/control), **Novel** (not in literature per lit-check), **Promising/Relevant** (disease/cell-biology relevance, testable). Priority = C×N×P. FINAL rule (aligned with REPORT.md): Tier 1 discovery ([!]) = audited C≥3.5 + independent replication + literature check (only K1 = KPNA4 S60); Tier 2 = audited/qualified C≈3 with cross-cohort/lab replication (K2–K4).

## Tree
- [~] 0. Infrastructure
  - [x] 0.1 Literature review → LITERATURE.md (6.9k words, ~120 refs), docs/lit/novelty_checklist.md (97 rows), docs/lit/SUMMARY.md
  - [x] 0.2 Compendium → results/infra_compendium/ (SUMMARY.md; loader scripts/infra_compendium/load.py; env /home/plm_env/bin/python)
    - 189,962 sites × 1,937 conds (4.6% filled); 850 inhibitor conds (522 high-conf targets), 172 time series, 86 dose series (decryptM = 790 conds).
    - Activity: core-v3 (503 kinases; evaluated curated_union library) + psp_z (189 kinases). Sanity: inhibitor targets decreased, median AUROC 0.95 (core-v3; optimistic: z_sig fitted in-sample) / 0.85 z_core / 0.86 psp_z.
    - Rules: group/CV by study_pmid; use z_core or psp_z for independent tests; filter n_substrates≥5; watch reversed / inhibitor-vs-inhibitor / background contrasts.
    - Audit-needed judgement calls: unnamed-compound targets (p38i, MEKi...), MRT ambiguity, PKD1kd→PRKD1.
  - [x] 0.3 CPTAC: harmonized phospho/protein/clinical/mutation, per-sample kinase activity → results/infra_cptac/ (SUMMARY.md)
    - 10 cohorts, 1,030 tumours + 523 NATs (492 paired); 99.7% site mapping, residue-verified; RNA/CNV/mut fetched from PDC; subtypes missing for CCRCC/HNSCC/OV.
    - Activity: core-v3 + mean-z on 4 inputs (tumour-centred, paired, protein-adjusted, paired-adjusted). Methods agree ρ≈0.53.
    - Sanity: MAPK1/3~RAS/BRAF, EGFR~amp (GBM .83), AKT~PIK3CA/PTEN, CDK~proliferation pass. ERBB2~amp FAILS (signal lands on ERBB3; ≤7 substrates) → never claim ERBB2 from activity; use pY sites/protein.
    - Caveats: z not calibrated (15–19% |z|>1.96, rank only); protein adjustment removes autophosphorylation signal; edge library differs from compendium's (to reconcile).
- [x] 1. Idea generation & ranking — docs/IDEAS.md (re-ranked after lit; wave-2 plan D2)
- [x] B-I08 Abundance vs PTM control in tumours → results/b_I08/SUMMARY.md
  - [x] I08.1 Activity decoupled from abundance; activating site beats protein (81/117 kinases) — partly known (Arshad 2019) C4/N3/P4
  - [x] I08.2 Genomics (mut/CNV) explains ~0 of decoupled kinase activity; other kinases' activity explains most (held-out ρ .78) — negative/useful
  - [x] I08.3 New activity-reporter sites: ABL1 Y185, CDK13 T588, PTK2 S840 (C3/N4/P4) — candidate; could check in compendium (do they move with inhibitors of their kinase?)
- [x] B-I17 Dose pharmacology → results/b_I17/SUMMARY.md (revised)
  - [x] I17.1 Downstream phosphosite potency shift fits ATP competition: implied cellular [ATP] 1.5–3.6 mM from 5 ATP-site drugs (Davis Kd + Patricelli Km); allosteric MEKi overshoot (replicated in independent KRAS dataset). C3–4/N3/P4 but n=5 drugs. decryptM 2.0 (133 drugs) restricted → asked user for access
  - [x] I17.2 Motif-predicted (not curated) target substrates are more potent (+0.33 log) C3/N3
  - [x] I17.4 Dasatinib→p38α/MK2 at ~1 µM (4 series; weak shuffle separation) C3/N2
  - [-] I18 cross-drug off-target assignment: ≈shuffle; scooped by decryptM 2.0 (Bayer 2025). I21 reach: restates "central kinases" — negative
  - [-] Ultrasensitivity vs network depth: not supported (down-sites steeper than up-sites, p=.005, fits I19 theme)
- [x] B-I19 Feedback/compensation atlas → results/b_I19/SUMMARY.md
  - [x] I19.1 (CANDIDATE key; [!] only after final audit — retro M1) **Primary-direction responses are portable, reverse-direction responses are cell-context-specific.** AUDITED (results/audit_I19/AUDIT.md): survives with qualification. Mirror test: in 52 stimulation studies increases replicate (log-odds 0.56→1.11→1.62 by effect size) and decreases don't (0.48→0.66→0.59; Δ −0.70 [−1.21,−0.25]); inhibitors show the reverse (Δ 0.30 [0.05,0.69] at |lfc| 1–2) but lean on MEK. Not MNAR, not noise (same-cell reverse changes replicate as well), gap opens when cell line changes. Original null had a bug (fixed). Wording: see audit. Auditor score C3/N4/P3 (inhibitor arm MEK-dependent; same-cell comparisons partly shared-control batch; 'two layers'/'compensatory' are interpretation). Fig: results/audit_I19/fig_mirror_inhibitor_vs_stimulation.png
  - [x] I19.2 Self-hyperphosphorylation: PAK4 S148 under PF-3758309 (3/3 cell lines, 1 study; novel?) ; AKT/PKC/MEK/ERK/CDK1 known
  - [x] I19.3 MET/ALK inhibition → Hippo/NDR (LATS1/2, STK38/L) activation, 6–7 studies (C3/N3) — candidate for follow-up (CPTAC: MET/ALK-high tumours vs YAP S127?)
  - [-] Most other kinase-level "up" hits: artefacts of 1–3 predicted substrates; classic transcriptional feedback not visible (≤1 h data)
- [-] B-I20 cis-CNV IV kinase→substrate in CPTAC → results/b_I20/SUMMARY.md — NEGATIVE: IV ≈ naive ≈ null except EGFR (pY adaptors SHC1, PTPN11 Y62, ERRFI1); converges with I08: kinase abundance rarely limiting in tumours (known, Arshad 2019). C2/N3/P2
- [x] B-I10/I06 Phospho-CMap + ICA modules → CPTAC → results/b_I10/SUMMARY.md
  - [x] I10.1 Consensus footprints transfer across studies (median AUROC .71); **down-site footprints transfer better than up-site (.68 vs .54) and only down-sites predict PIK3CA/PTEN in tumours → independent replication of Finding A, incl. tumour side**
  - [x] I10.2 CANDIDATE B — AUDITED, QUALIFIED (results/audit_I10/AUDIT.md): footprints coherent in 10/10 cohorts, robust to protein adj, co-abundant/complex pair removal, strict matched nulls (z≈20–30), top-5-PC removal (MEKi .153→.068), ≈4× typical PSP set. BUT −0.88 MEK/EGF partly built in (sign composition); axis = MAPK1/3 output; **~half shared with a cold-ISCHAEMIA footprint (|r| .5–.75 all cohorts) → possibly pre-analytical**; DNA-damage footprint cleanest (unaffected by ischaemia/PCs). Asymmetry confirmed as PRIMARY-direction (inhibitor down 75%, stimulus up 75%). Auditor: qualified.
  - [-] I10.3 Immune ↔ PTP-inhibitor mimicry FAILS as stated (CPTAC barely measures pY; driven by pS/pT HSPB1 S15, AHNAK, AKAP13, MARK2 + STAT3 Y705; not PTP-specific)
  - [-] Genotype recovery no better than kinase scores; ICA: no novel module; survival null; IFN footprint fails
- [x] A-I19 Independent audit of up/down asymmetry (MNAR, base rates, time/dose, study composition, mirror test in stimulations) → results/audit_I19/ (4 thr)
- [x] B-I05 Timing/directness/portability → results/b_I05/SUMMARY.md
  - [-] Reverse layer is NOT later (inhibitors: decreases accumulate; timing doesn't explain I19 asymmetry)
  - [x] Portability tracks persistence, not onset
  - [x] Supports I19.1: motif×direction interaction (stimuli, Δ .94 [.41,1.42], beats random sets); reverse responses ~2× less potent (decryptM, 18 drugs); reverse kinase signature (NEK, IRAK4, MAP4K/STE20, MAP3K, MARK/SIK) cell-line dependent
- [x] V1 Cross-validation → results/val_1/SUMMARY.md
  - [-] I19.3 MET/ALK→Hippo: NOT supported in CPTAC (RTK/MAPK axis positively assoc. with Hippo kinases, r +0.19–0.23)
  - [-] I19.2 PAK4 S148: unreplicated (no other PAK4 inhibitor); likely nonspecific (also up in mitosis, ulixertinib, abemaciclib, vemurafenib 24h); activation-loop S474 goes down
  - [-] I08.3 reporter sites: only IKBKB S697 validated (IKKi AUROC .81, TNF .75); ABL1 Y185 not supported; CDK13/PAK2/PTK2 unreplicated (untestable) or family-level
  - NOTE: stored activity scores include autophosphorylation edges → use masked rescoring (results/val_1/tables/corev3_masked_long.parquet) for any site-vs-own-kinase analysis
- [-] B-I01 Dark kinases → results/b_I01/SUMMARY.md — NEGATIVE for biology: 66% of dark-kinase activity is family echo; family-disjoint residuals pass drug-target control for well-studied (AUROC .60) but not dark kinases; no context/partner survives FDR or beats family-label baseline. Method caution only.
- [x] B-I22 Is reverse layer cell-intrinsic? → results/b_I22/SUMMARY.md — NO: reproducible within a study's culture of a cell line (Hijazi), not across labs for the same cell line (MCF7 Hijazi vs Ochoa); baseline RNA/proteome/DepMap don't predict it (positive control passes). Reframes Finding A as a cross-lab REPRODUCIBILITY asymmetry (primary portable; reverse study×culture-specific).
- [-] B-I23 Reverse layer vs drug resistance → results/b_I23/SUMMARY.md — NEGATIVE (pre-registered; opposite sign, fragile, no CV value; underpowered: 3 lines, 60 min)
- [-] B-I13 Phosphatases → results/b_I13/SUMMARY.md — NEGATIVE: PPP2R2A loss, PPP2R1A hotspots, PPM1D gain, PTP loss don't shift rule-matched sites (positive control PTEN→PRAS40 T246 10/10 passes); phosphatase lesions don't bias kinase scores; compendium too thin for rule learning (OA vs cantharidic r .21). C2/N4/P2
- [x] A-I10 Audit candidate B (footprint tumour coherence: protein/PC/purity confounds, matched nulls) + C (immune–PTP mimicry vs leukocyte pY content) → results/audit_I10/ (4 thr)
- [x] B-I03 New kinase→site edges → results/b_I03/SUMMARY.md
  - [-] Global ranking: new evidence does not beat kinaseSub+PPI+motif
  - [x] Confirmation layer: replicated decrease under kinase LoF raises P(edge) ~2.5× among high-prior pairs (TEST proteins OR 2.3, external OR 2.8; shuffles null); family-level only
  - [x] Rediscovered literature-only MTOR→AMPKα1 S347 (positive control); conflicts: CK2 (not MAPK)→KPNA4 S60 (acidic motif), AKT→ARHGEF12 S1288
  - [x] CANDIDATE D (validated by V-I03 → results/val_I03/SUMMARY.md, tables/verdicts.tsv):
    - LAMTOR1 S98 likely DIRECT mTORC1 (catalytic mTORi z −5.1; S6Ki ratio ≈0 in 3 lines; pEC50 matches PRAS40 S183) — novel site/kinase; against: AKTi-sensitive, no tumour co-variation
    - RRAGC S15 probably direct mTORC1, provisional (rapalog/catalytic ratio .48 like 4E-BP1 S65; PI3Ki no effect while S6K site fell 4.5×) — next to known S21 feedback site
    - **CK2→KPNA4 S60** (CX-4945 z −2.5…−7.0 in 3 lines; 0/25 ERK-inhibited exps; acidic motif) — corrects Oncogene 2019 "MAPK-dependent"
    - ERK pathway→EPN1 S435 (r .57 with ERK A-loop; MEK/RAFi −4…−10; dose matches ERK) — direct unproven
    - [-] TBC1D5 S730: not mTORC1-specific (AGC/MAPK input)
- [-] B-I27 Allele-specific signalling → results/b_I27/SUMMARY.md — pre-registered P1–P5 not confirmed/underpowered (P1 G12R AKT deficit: right sign, p .25, MDE .67 SD; PIK3CA helical vs kinase unanswerable in CPTAC). Exploratory: EGFR-mut LUAD ~1 SD lower ERK output at similar receptor pY (PI caveat: WT group contains KRAS-mut — needs RAS/RAF-WT comparator; possibly supports E); PDAC G12R MAPK excess vs G12D (q .015); BRAF V600E COAD looks RTK/PI3K-inhibited (q .015, supports P4 via footprints). C3/N3/P3
- [x] B-I26 Perturbation audit of PSP annotations → results/b_I26/SUMMARY.md (C3/N3/P3)
  - Only 680/9,305 edges testable (7%); 17% supported vs 5.3% random (3× enrichment); well-established edges 22% supported → 'unsupported' ≠ wrong
  - CANDIDATE F: 27 curated corrections; dominant pattern AKT/CDK/PKC-annotated sites behaving as ERK/RSK sites; literature-confirmed: CHEK1 S280 RSK not AKT, EPHA2 S897 RSK, CTTN S405 ERK; new ERK/RSK candidates DBNL S269, TFPT S180, AKAP11 T1100
  - Pruning unsupported edges: +0.005 target AUROC (p .005) — small
  - !! CONFLICT: CX-4945 does not lower known CK2 sites at ~1h → V-I03 'CK2→KPNA4 S60' may be CX-4945 off-target (CLK/DYRK/PIM/HIPK) — resolve in final audit of D
- [x] B-I24 Cascade gain in tumours → results/b_I24/SUMMARY.md
  - [x] I24.1 Lesions shift OFFSET not GAIN: 0/162 gain tests (power 84–94% for ~35–40% change); KRAS/PTEN/TSC offsets as expected — powered negative, C4/N3/P3
  - [x] I24.2 CANDIDATE E (AUDITED, reworded; results/audit_I24/AUDIT.md): "non-transitivity" FAILS as worded (only unadjusted data; vanishes with RTK protein covariate). Defensible: **abundance-driven receptor pY does not reach ERK in tumours — at fixed EGFR/ERBB2 pY, more receptor protein → lower ERK (−0.15 SD [−.25,−.05], p .002, 6/7 cohorts); acute receptor perturbations in cell lines are transitive (RTK→ERK .58 vs expected .51)** → chronic abundance buffered, acute activity propagates. Noise simulation predicts opposite (P .01). Auditor C3/N3/P4. Fig: results/audit_I24/figures/AF1_decomposition_specifications.png
- [x] A-I24 Audit done
- [-] B-I25 Mechanism of Finding A → results/b_I25/SUMMARY.md — NEGATIVE (all 4 pre-registered H fail): reverse layer not cell-cycle, energy/stress, Hippo, serum, normalisation; drug-induced mitotic shifts are PORTABLE (class-explained). Lab-specific reverse component has no identified mechanism (shared within-run structure untested) → Finding A stays a reproducibility statement. STOP investing in A (lesson 11).
- [x] B-Q8 Tumour vs NAT wiring → results/b_Q8/SUMMARY.md — strengthens B: cell-line perturbation wiring is preserved in NORMAL tissue and tumours (MAPK axis = PC1 in every NAT cohort; 529 signed compendium relationships 91% sign agreement in both); only robust tumour change = GAINED coherence of DNA-damage/replication-stress programme (meta z 5.3–6.3, survives composition/ischaemia/proliferation); apparent RTK loss = PC artefact. C4/N3/P3
- [x] A-D Audit of candidate D → results/audit_D/AUDIT.md (auditor C/N/P):
  - D3 **KPNA4 S60 = fast-turnover CK2-type site, not MAPK** — SURVIVES w/ qualification (3.5/3/3.5): CK2 inhibition verified in Rusin 2017 mitotic HeLa (41% known CK2 sites fall; S60 rank 63/12,455); member of a 14-site fast CK2-responsive set (64% strong CK2 motif, OR 16, p 1.5e-6); CX-4945 off-target motifs excluded; 0/25 ERK-inhibited exps. Resolves I26 conflict (Hijazi 1 h: CK2 set doesn't fall). Needs SGC-CK2-1/rescue.
  - D1 LAMTOR1 S97/S98 (localisation ambiguous) candidate direct mTORC1 — QUALIFIED (3/4/3): z −5.07 (0.11th pct), LOSO −4.9…−6.1, other sites flat; but both studies Cutillas lab, S6Ki not verified, AKTi ratio .35
  - D2 RRAGC S15 — QUALIFIED/weaker (2.5/4/3): mTOR motif 100th pct, 6/6 organoid lines; ~1.5 usable studies
  - D4 EPN1 S435 ERK-pathway-dependent (3/3/3); direct untested
- [x] B-I28 Turnover atlas → results/b_I28/SUMMARY.md — CANDIDATE H: turnover is a reproducible site property (split-half .34 vs .16), poorly predicted by sequence/structure (AUROC .545 vs .725 same-site), phosphatase-consistent determinants (pThr faster, acidic stable); fast sites are tumour reporters (.27 vs .13); fast-edge library +.032 AUROC on edited held-out studies. Self C3.5/N4/P3.5 → audit A-H
- [x] A-H Audit → results/audit_H/AUDIT.md (auditor C3/N2.5/P2.5): H1 qualified (lab split .27 vs .15; specific to own-kinase inhibitors; ERK/MTOR carry it, CDK fails); H2 pThr real within proteins but CDK/mitosis-driven & known (B55); acidic/basic determinants FAIL (same in non-substrates); H3 gap halves after (over-)adjustment (.117→.063); H4 survives leave-one-lab-out (+.035) but carried by few studies; condition-weighted +.015
- [-] B-I29 CANDIDATE G (ERK ceiling in EGFR-mut LUAD) → results/b_I29/SUMMARY.md — robust in CPTAC LUAD discovery (measured ischaemia, country, batch, 2 pipelines; GBM mirrors) but FAILS replication (CPTAC LUAD Confirmatory −0.03; Taiwan +0.17 underpowered). Cohort-specific observation only. Incidental: compendium ischaemia footprint does not track measured GDC ischaemia time in LUAD → B's ischaemia caveat to be re-tested with real times (V-B).
- [x] V-E Ischaemia re-check → results/val_E/SUMMARY.md: E abundance coef −.151→−.158 (p .001, 6/7), KRAS offset +.67→+.60, LUAD EGFR-mut ERK −.76→−.77 — all ROBUST; ischaemia score uncorrelated with receptor protein/KRAS (|ρ|≤.14). Scores saved: results/val_E/ischaemia_scores.tsv
- [x] A-F Audit → results/audit_F/AUDIT.md
  - [!] **KPNA4 S60 (importin-α3) is CK2α-dependent — public chemical-genetic confirmation** (Menyhart 2023, PXD038050, U2OS): SGC-CK2-1 4 h −3.73 log2 (FDR 3e-4), rescued +2.98 (~80%, FDR 4e-4) by inhibitor-resistant CK2α; CX-4945 −3.1. All 6 measured members of audit_D's fast CK2 set fall and are rescued. Corrects Hazawa 2019 'MAPK-dependent'. Open: direct vs indirect (acidic motif supports direct), 1-h kinetics, other cell types.
  - [x] F downgraded (auditor C2/N2/P2.5): AKT-annotated basophilic sites depending on MEK/ERK beat chance (10/30 vs 7.6%, p .002) but RSK vs AKT not separable (AKTi relieves feedback onto ERK/RSK; motifs ρ .82); CDK half at chance (only 1 lab verified). EPHA2 S897/CHEK1 S280 already PSP-RSK; MAP2K2 T394 contradicted (Cdk5, Banks 2015); new candidates DBNL S269, AKAP11 T1100, PIK3C2A S259. Labs ≠ studies (Cutillas, Kuster).
- [x] V-B Measured ischaemia → results/val_B/SUMMARY.md: GDC excision-to-freezing times (7 cohorts, n=723) uncorrelated with MEK/EGF/TNF footprints (pooled |ρ|≤.09) and with the compendium 'ischaemia score' (ρ −.03); adjusting for measured time leaves coherence unchanged (MEK .170→.169). The ischaemia footprint is a co-readout of the MAPK/stress axis → audit_I10's ischaemia caveat WITHDRAWN (over-adjustment). B restored (residual: 3 cohorts without times; warm ischaemia untested).
