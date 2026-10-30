# I23 pre-registration (written 2026-10-29 before any drug-sensitivity data was downloaded or viewed)

Question: does the reverse-layer (increase) phosphoproteome response of a drug x cell pair predict resistance, beyond cell identity/baseline and primary-inhibition magnitude?

## Primary test (P1) — one test, declared now
- Data: Hijazi 2020 (PMID 31959955), single study, 60 min, same dose per drug in all 3 cell lines (HL60, MCF7, NTERA2) -> dose/time/study controlled by design. Single-agent, no background agent.
- Preprocessing: log2FC per condition; within study x cell, subtract each site's median across all inhibitor conditions of that cell (I19 centring). Sites measured in the condition only.
- Exposure: reverse/primary log ratio LR = ln((f_up + 0.002)/(f_down + 0.002)), f_up/f_down = fraction of measured sites with centred log2FC >= +1 / <= -1.
- Outcome: GDSC2 AUC (fallback order if a drug x cell is missing: GDSC1 AUC, then PRISM secondary AUC; each source z-scored within drug before pooling is NOT done for P1 — P1 uses the single best source with most matched pairs, chosen by match count only, not by results). Higher AUC = more resistant.
- Model: AUC ~ drug FE + cell FE + LR (OLS). Drug FE removes drug potency; cell FE removes global cell sensitivity/baseline/lineage. The LR coefficient tests drug x cell-specific association.
- Hypothesis: beta(LR) > 0. Inference: two-sided permutation p (LR permuted among the 3 cells within each drug, 10,000 perms) and drug-cluster bootstrap 95% CI. Report standardized beta (per SD of residualized LR).
- Requires >= 15 drugs with all 3... (>= 2 cells matched per drug); else declared underpowered and reported descriptively.

## Secondary (declared, exploratory-labelled)
S1. Same model with log f_up and log f_down as separate terms (reverse controlling for primary magnitude).
S2. Add target-engagement covariate (mean core-v3 z_core of the drug's curated primary targets in the condition) and baseline target expression (DepMap, if obtained).
S3. Pre-specified kinase groups in the reverse layer (core-v3 z_core, positive part, per condition): STE20/MAP4K (MAP4K4,TNIK,MINK1,MAP4K1-5), NEK (NEK1-11), IRAK4, MAP3K (MAP3K*), MARK/SIK (MARK1-4,SIK1-3), AKT (AKT1-3), SFK (SRC,YES1,FYN,LYN,LCK,HCK,FGR,BLK). Each: AUC ~ drug + cell + z_group; BH across 7; leave-one-drug-out sign stability.
S4. Other datasets (decryptM multi-cell drugs, crizotinib/lorlatinib neuroblastoma, vemurafenib COLO205/HCT116): within-drug descriptive correlations only.
S5. Negative/controls: same model with LR from a random 50% site subset; shuffled exposure; primary magnitude alone (log f_down).
