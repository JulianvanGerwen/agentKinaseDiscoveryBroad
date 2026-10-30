# IDEAS — pool and ranking

Scores 1–5. **Feas** = feasible in ≤1 day with local data; **Nov** = expected novelty (pre lit-check); **Imp** = impact/relevance; **Disc** = chance of producing a compelling positive result. Priority = Nov×Imp×Disc (Feas gates scheduling).
Status mirrors RESEARCH_TREE.md. Pre-literature draft (T+0); re-ranked after wave 1.

| ID | Idea | Data | Feas | Nov | Imp | Disc | Pri |
|---|---|---|---|---|---|---|---|
| I01 | Dark-kinase function by guilt-by-association: co-variation of activity (and of their own regulatory phosphosites) across 1.9k conditions + CPTAC → assign understudied kinases to pathways; validate with interactomes (BioPlex/OpenCell), DepMap co-dependency | comp+CPTAC+DepMap | 4 | 4 | 4 | 3 | 48 |
| I02 | Systematic feedback/rewiring atlas: for each inhibited kinase, which other kinases go up (compensation) consistently across studies/cell types; distinguish known (MEK→RTK, mTOR→AKT) from novel | comp (inhib) | 4 | 3 | 5 | 4 | 60 |
| I03 | Orphan sites → kinase by co-regulation + motif: sites with no known kinase whose profiles track a kinase's activity across conditions AND match its Cantley motif → new kinase–substrate assignments | comp+Cantley | 4 | 3 | 3 | 4 | 36 |
| I04 | Activity-reporting sites: which kinase autophosphorylation/regulatory sites truly track inferred activity (compendium+CPTAC); find new "activity marker" sites | comp+CPTAC+PSP reg | 5 | 3 | 3 | 4 | 36 |
| I05 | Temporal ordering / lead–lag causality in time series (Granger / transfer entropy borrowed from econometrics) → directed kinase–kinase network; validate with inhibitor perturbations | comp (time) | 3 | 4 | 4 | 2 | 32 |
| I06 | Signaling modules by ICA/NMF of the site×condition compendium (as in transcriptomics) → modules not explained by known kinases = unknown regulators; module activity in CPTAC vs survival | comp+CPTAC | 4 | 3 | 4 | 3 | 36 |
| I07 | Pan-cancer kinase activity vs somatic mutations beyond canonical drivers (e.g., which non-kinase mutations rewire kinase activity); protein-adjusted | CPTAC | 4 | 2 | 5 | 3 | 30 |
| I08 | Post-translational vs expression control: for each kinase in tumours, is activity explained by its protein abundance, mRNA, or neither → "PTM-controlled" kinases | CPTAC | 5 | 4 | 3 | 4 | 48 |
| I09 | Kinase activity predicts survival / therapy-relevant states consistently across cohorts (meta-analysis), protein-adjusted, vs proliferation confound | CPTAC | 4 | 2 | 5 | 3 | 30 |
| I10 | Perturbation-signature transfer: map tumour phospho profiles onto compendium perturbation signatures (connectivity-map style, borrowed from CMap) → which drug/stimulus "states" tumours resemble; candidate drug repurposing | comp+CPTAC | 4 | 4 | 5 | 3 | 60 |
| I11 | Off-target/indirect effects: inhibitor-driven kinase activity changes not explained by Klaeger binding → network-propagated vs novel off-targets | comp+Klaeger | 4 | 3 | 4 | 3 | 36 |
| I12 | Context-specificity: how much does a kinase's substrate response depend on cell type (substrate "rewiring" across contexts) — variance partitioning of site responses to the same perturbation | comp | 3 | 4 | 3 | 3 | 36 |
| I13 | Phosphatase activity inference (reverse logic) and kinase–phosphatase balance across conditions/tumours | comp+CPTAC | 3 | 4 | 4 | 2 | 32 |
| I14 | Evolutionary/functional-score view: are sites that respond to many perturbations more functional/conserved (Ochoa 2020)? "Hub" vs "private" sites | comp+Ochoa | 5 | 2 | 3 | 4 | 24 |
| I15 | Kinase activity vs DepMap dependency: tumour kinase activity outliers whose cell-line dependency is selective → activity-based target nomination | CPTAC+DepMap | 3 | 3 | 5 | 2 | 30 |
| I16 | Sex/age differences in kinase activity across cancers | CPTAC | 5 | 3 | 3 | 2 | 18 |

## Re-ranking after literature review (T+~1h)
Novelty downgrades: I01 (co-regulation→function is known; keep only DepMap/interactome-validated dark-kinase angle, Nov 3→2), I07/I09/I15 (CPTAC driver/survival/DepMap claims largely taken, Nov 2), I14 (circular with Ochoa functional score, drop).
New ideas from literature gaps:

| ID | Idea | Data | Feas | Nov | Imp | Disc | Pri |
|---|---|---|---|---|---|---|---|
| I17 | Dose-pharmacology map: per-site EC50 (decryptM dose series) vs target Kd_app (Kinobeads) → amplified/ultrasensitive sites; does network distance from target predict EC50 shift & Hill slope? | comp (decryptM)+Klaeger | 4 | 5 | 4 | 4 | 80 |
| I18 | Functional off-targets: site EC50s that track a non-intended kinase's Kd_app across drugs | comp+Klaeger | 4 | 4 | 4 | 3 | 48 |
| I19 | Self-hyperphosphorylation & feedback invariance atlas: kinase regulatory sites and other sites that RISE under inhibition, conserved across cell lines/drugs (merges I02) | comp | 5 | 4 | 4 | 4 | 64 |
| I20 | Mendelian-randomisation-style causal kinase→substrate tests in CPTAC (cis-CNV/truncation as instrument for kinase protein) → validates/discovers edges in vivo | CPTAC | 4 | 5 | 4 | 3 | 60 |
| I21 | Kinase "reach": phosphoproteome fraction perturbed per unit target occupancy vs centrality/essentiality | comp+Klaeger+DepMap | 3 | 5 | 3 | 3 | 45 |
| I10 | Phospho-CMap (kept): tumours scored against proteome-scale perturbation footprints | comp+CPTAC | 4 | 4 | 5 | 3 | 60 |
| I06 | ICA phosphomodules on compendium → projected onto CPTAC (kept, Nov 3→4 with CPTAC projection) | comp+CPTAC | 4 | 4 | 4 | 3 | 48 |

**Wave-2 plan (decision D2):** launch I17+I18+I21 (dose pharmacology), I19 (feedback/self-hyperphosphorylation), I10+I06 (CMap/modules→tumours), I20 (MR in CPTAC) as compute frees. I08 already running. Rationale: highest Nov×Imp, all feasible with local data, and they span different data types (dose series, inhibitor panels, tumours) so failures are uncorrelated.
