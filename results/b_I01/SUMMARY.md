# I01 SUMMARY — dark kinases net of family signal (saved by PI from agent message)
Design: kinase-specific substrate sets (top-1% for k, not top-1% for family + 5 motif neighbours), residual score minus confounder score, strict pc10 variant (PC1 of kinome activity = 69%); split-half coherence vs same-group random sets; drug-target positive control; orthogonal validation (DepMap, interactions, GTEx, localisation) vs family-label baseline.
1. **Family echo quantified:** family+neighbours explain median R² .66 of a dark kinase's all-substrate activity (.58 well-studied); specific score .20.
2. Specific signal coherent (split-half .114 vs .034 random, p 1e-15; 26/123 dark pass FDR; pc10: VRK2, TLK1, MAPK15, DSTYK, DYRK4, PRP4K, ADCK5).
3. **Drug control:** residual falls when drug binds k but not family — well-studied AUROC .596 (n 80, p 4e-6), family score unmoved; dark .52 (n 21) not validated.
4. Residual still echoes motif class (DYRK3~RSK/AKT, TLK1~PLK).
5. No dark-kinase context survives FDR (min .26; positive controls only nominal).
6. Orthogonal validation weak: DepMap co-dependency AUROC .567 for 26 dark kinases but vanishes for variably-essential kinases and is beaten by family-label baseline (.59–.82); same-band CNV artefacts. Kinoplex files are a phosphosite atlas, not interactome.
Candidates (unsupported): VRK2 falls under CDK inhibitors (likely direct drug binding); TLK1–PLK1 (motif echo).
Novelty: CoPheeKSA 2025 annotates understudied kinases; family-echo fraction + drug control new. Self C2/N3/P2.
Files: scripts/b_I01/01–09; results/b_I01/tables, figures fig1–3.
