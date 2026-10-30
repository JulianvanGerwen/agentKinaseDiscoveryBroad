# Literature summary for idea generation (≤1 page)
(Written by the literature agent; saved by the PI. Full review: /workspace/LITERATURE.md; known findings: docs/lit/novelty_checklist.md, 97 rows.)

## 10 takeaways
1. The bottleneck is annotation, not algorithm: <5% of sites have a known kinase, ~90% of annotations on 20% of kinases, >150 kinases without substrates (Needham 2019, PMID 30670635); library > statistic (benchmarKIN, PMID 40404650). Annotation-free approaches (perturbation footprints, motifs, modules) have most room.
2. Kinase-level compendium atlases exist: Ochoa 2016 (215 kinases × 399 conditions, PMID 27909043), Invergo 2020 signed kinase–kinase network (PMID 32437683), Hijazi 2020 inhibitor topology (PMID 31959955). "Co-regulation predicts complexes/substrates" and "central kinases regulated most" are NOT novel. Site-level modules and dose structure are open.
3. Most phosphosites weakly constrained; ancient sites at interfaces are functional. Ochoa 2020 functional score uses regulation as a feature → "functional sites are more regulated" is circular.
4. CPTAC kinase claims taken: driver cis-effects; CDK4/6–Rb and proliferation CDKs; DNA-damage kinases in HRD (Geffen 2023); PTPN11 downstream of RTKs; immune-subtype kinases (Petralia 2024); outlier kinases vs DepMap (Gillette 2020, Savage 2024); drivers converge on shared kinase states (Li 2023 Cell, PMID 37582357).
5. CPTAC co-phosphorylation is dominated by protein abundance and purity → any tumour claim needs abundance/purity adjustment and per-cohort replication.
6. Known feedback is axis-specific (MEK→MYC→RTK, AKT→FOXO→RTK, RAF paradox, mTOR→AKT). No dose-resolved feedback atlas across contexts.
7. decryptM (790 of our contrasts) + Kinobeads Kd_app under-used: site EC50 vs target Kd separates direct, amplified and off-target effects.
8. Data-only network inference is weak (HPN-DREAM); claims need interventional data and prior-knowledge baselines.
9. "Cell-line kinase activity predicts drug sensitivity" is published (Koenig/Olsen 2026 NSMB, PMID 42816593; Frejno 2020) — avoid as headline.
10. Motif-based activity is family-level; inhibitors can raise their own target's activation-loop phosphorylation (Okuzumi 2009).

## 15 open questions
1. Dose-sensitivity map: site EC50 vs target Kd_app (decryptM × Klaeger); amplifier sites; network depth vs EC50 shift / Hill slope (ultrasensitivity proteome-wide untested). *Verify decryptM supplement before claiming novelty.*
2. Functional off-targets: site EC50 tracks a non-intended kinase's Kd_app.
3. Feedback invariance atlas: sites rising on inhibition, conservation across 13 decryptM + 3 Hijazi cell lines.
4. Catalogue of self-hyperphosphorylation of kinases under their own inhibitors (known only AKT, RAF).
5. ICA phosphomodules on compendium projected onto CPTAC.
6. Phospho-CMap: tumours scored with proteome-scale inhibitor footprints (P100 CMap had ~96 sites).
7. MR-style kinase→substrate causal effects in CPTAC (cis-CNV/truncation as instrument).
8. Do perturbation-derived kinase–kinase edges hold in tumours? Rewiring cell lines vs patients.
9. Per-site variance partition in CPTAC (abundance vs kinase activity vs residual).
10. Dark-kinase activity atlas vs DepMap co-essentiality and GTEx.
11. Site evolutionary age vs dynamic behaviour (response breadth, dose sensitivity, tumour variance).
12. Kinase "reach": fraction of phosphoproteome perturbed per unit target occupancy vs centrality/essentiality.
13. Kinase hierarchy from footprint nesting (nested effects models) on inhibitor panels.
14. Co-regulated but not co-essential kinase pairs → redundancy/compensation.
15. Temporal ordering across studies: are early robust responders enriched for direct-kinase motifs in human data?

Borrowed ideas (mapped in LITERATURE.md §6): ICA/iModulons, MOFA, Mendelian randomisation, modular response analysis, CMap/PROGENy footprints, network control theory, DebtRank, ecological nestedness, pseudotime over dose, Pareto archetypes, nested effects models.

Caveats: Koenig/Martinez-Val/Olsen NSMB 2026 (33-cell-line atlas) is a novelty threat for cell-line activity–drug-sensitivity claims. The GitHub "awesome-pan-cancer-phosphoproteomics" list contains fabricated citations — do not use.
