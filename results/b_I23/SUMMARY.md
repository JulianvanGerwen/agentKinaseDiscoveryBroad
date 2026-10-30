# I23 SUMMARY — reverse layer vs drug resistance (saved by PI from agent message)
Pre-registered (PREREG.md). Hijazi 2020: 32 drugs × HL60/MCF7/NTERA2 (94 pairs, 60 min, same dose); GDSC2/1 AUC; LR = ln(frac up≥1 / frac down≤−1); AUC ~ drug FE + cell FE + LR; within-drug permutation, drug-cluster bootstrap.
- **Hypothesis refuted in sign:** LR→AUC std β −0.28 [−0.47,−0.05] (p .025; larger reverse layer = slightly MORE sensitive), fragile: n.s. dropping HL60 or MCF7, or 3 drugs (MEKi in HL60, amuvatinib); GDSC2-only p .10; 18/14 per-drug sign split.
- No cross-validated predictive value (LODO R² .03); no kinase group passes FDR (STE20/MAP4K −0.20 p .12; AKT +0.17 p .13).
- Underpowered: ~8 informative drugs, 3 lines; 60-min data precede adaptive feedback (Duncan 2012, Lito 2012).
Novelty: up vs down decomposition untested before, but result negative. Self C2/N3/P2.
Files: results/b_I23/PREREG.md, tables T0–T9, figures/fig1; scripts/b_I23/01–07.
