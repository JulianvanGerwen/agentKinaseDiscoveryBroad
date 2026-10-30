# I24 SUMMARY — cascade signal-transmission gain in tumours (saved by PI from agent message)
Setup: 10 CPTAC cohorts, 1,030 tumours; stage composites of direct sites, no shared sites: RTK pY → adaptor pY (SHC1, GAB1, PTPN11) → ERK1/2 → RSK / ERKSUB (479 PSP ERK substrates); adaptor → AKT substrates → mTORC1 (4EBP1) → RPS6. Huber slopes adj. purity + E2F, RE meta. MEK/AKT/S6K activation loops unmeasured.
1. Gains reproducible (all core edges positive 10/10 cohorts; standardized 0.21–0.41); noise attenuation ~1.5–2× (IV gains .55–.75).
2. **Non-transitivity (candidate E):** RTK→adaptor .41, adaptor→ERK .41, but RTK→ERK .10 [.02,.18] (noise-corrected .12). 2SLS decomposition: RTK-driven adaptor pY → ERK .12 [−.05,.30] vs RTK-independent .46 [.34,.57], p 4e-4 (7 cohorts, unadjusted); AKT substrates .22 vs .40 (p .04); protein-adjusted weaker (.28 vs .44, p .07). Generalises EGFR-amp GBM decoupling. Needs weak-instrument check. C3/N4/P4.
3. **Lesions shift offset, not gain (powered negative):** 0/162 pre-declared gain tests (DUSP4/5/6, SPRY, SPRED, NF1, RASA1, PTPN11, PTEN, INPP4B, PHLPP, TSC, STK11, PPP2R2A, KRAS/BRAF/EGFR/PIK3CA) pass (min q .99); spike-in power 84–94% for 35–40% gain change (DUSP6, PTEN); genome-wide scan at null rate (17q11.2 fragile). Offsets as expected: KRAS mut RTK→ERK +0.67 SD (q .005), PTEN loss, TSC1/2 protein; DUSP6 protein −0.19 (q .15). C4/N3/P3.
4. Controls: residuals not driven by downstream protein; unadjusted ≈ adjusted; permutation-calibrated.
5. LUAD NMF C2 (KRAS/STK11/KEAP1) highest adaptor→ERK residual (+0.93 SD, q .006); survival null.
Novelty: KRAS/PTEN/TSC offsets known; DUSP4→ERK known (Balko 2012); gain map, IV non-transitivity, offset-not-gain new. 
Files: scripts/b_I24/01–09; results/b_I24/T1–T8, F1–F3.
