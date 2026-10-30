# I03 SUMMARY — new kinase→site edges from integrated independent evidence (saved by PI from agent message)
Setup: 63 kinases with LoF (inhibitor ≥100 nM curated target or genetic) in ≥2 studies; 288,583 evaluable pairs. Evidence: replicated decreases (primary layer), dose coherence (9 kinases), tumour co-variation with activating site (34 kinases), PPI, Cantley motif, kinaseSub (release; and DEV no-phospho variant for held-out TEST proteins).
Validation: new evidence alone weak (TEST AUROC .57, external .59); does NOT improve global ranking over kinaseSub+PPI+motif (Δ −0.001). **As a confirmation layer on high-prior pairs (pct≥95): replicated decrease OR 2.3 [1.4,3.9] on TEST proteins, 2.8 [1.7,4.7] on external non-PSP edges**; tumour z>3 OR 2.0 (TEST); dose r>.3 OR 5.3 (external); shuffles ≈1.0. Family-level only (within-family label swap AUROC .577 vs .594; top-1 exact 40%, family 54%).
Nominations (675; family-level; filtered vs PSP/OmniPath/iPTMnet; explaining-away; no MEK/RAF):
- **MTOR→PRKAA1 S356 (AMPKα1 S347)**: rediscovered literature-only edge (Ling 2020 Nat Metab) — positive control.
- **mTORC1 lysosomal cluster**: MTOR→RRAGC S15 (3/3 studies), LAMTOR1 S98 (2/2; tumour 3.2), TBC1D5 S730 (2/3; tumour 3.9). Novel sites (RagC S21 known mTORC1 feedback, Yang 2019).
- ERK1/2→EPN1 S435 (3/3; tumour 5.7; MEKi 9/10 studies); SAMD4B S271; LARP1 T376.
- AKT→ARHGEF12 S1288 (reported RSK2 site — conflict/extension); AKT→ZNF185 S64.
- CK2→KPNA4 S60 (paper says MAPK-dependent; MEKi no effect here — conflict); CK2→FXR2 S453.
- ATM→TP53BP1 S580; ATR→MCM6 S13, CLSPN S950, ACIN1 S243; ROCK→PPP1R12A S903; others.
Self C3/N3/P3. Caveats: downstream-kinase confound (MTOR via S6K/AKT), CDK1 mitotic-exit confound, dose only decryptM.
Files: results/b_I03/tables (top20_curated_edges.tsv, edges_novel_candidates.tsv, edges_ranked_all.tsv.gz, val_*), figures fig1–3; scripts/b_I03/01–09.
