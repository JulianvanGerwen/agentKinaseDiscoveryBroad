# I20 SUMMARY — cis-CNV instrumental-variable tests of kinase→substrate edges (saved by PI from agent message)
Design: 10 cohorts, 1,013 tumours; instrument = kinase CNV (resid. purity, CIN); outcome = protein-adjusted site; RF + Wald/2SLS, IV meta where F>10; controls: 400 non-kinase fake instruments, permutations, focal instrument, co-amplified neighbour adjustment.
1. Instruments strong (F>10 in 993/3,191 kinase×cohort; 256 kinases) but kinase protein is a weak exposure.
2. **Negative: IV ≈ naive correlation ≈ barely above null** (PSP AUROC .525 vs naive .535, null .501; motif .502); signal mostly arm-level co-amplification (vanishes with focal/neighbour adjustment); 8.7% of kinases beat fake-instrument null (5% expected).
3. EGFR only clean positive: pY-specific (12/20 top sites pY vs 1.2% bg): PLCG1 Y783, ANXA1 Y21, GAB1 Y627 (PSP) + SHC1 Y349/Y427, ERRFI1 Y394, PTPN11 Y62 (z 9.3; SFK-mediated per 2025 preprint → indirect). Robust to focal/neighbour adjustment.
4. No confident new edges: motif-step FDR ≈ 0.6 (18.7% vs 11.5% swapped-motif null). Candidates EGFR→CFL1 Y140, TLN1 Y70, CSTB Y97, AHNAK Y715.
5. PSP edges: 76 supported, 56 contradicted (≈34 each by chance) of 1,378 → "kinase abundance not limiting", not wrong annotations.
Novelty: closest Memon 2021 Cell Rep Med (RPPA CNV→phosphosite); IV-vs-naive comparison new but negative. Self C2/N3/P2.
Files: scripts/b_I20/00–06; results/b_I20/tables T1–T6, figures fig1–4; cache ~1.5 GB deletable.
