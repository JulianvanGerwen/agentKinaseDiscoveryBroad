# b_I22 notes — pre-declarations (written 2026-10-29 before any b_I22 analysis was run)

Prior knowledge before declaration: I05 T5c (read before declaring) shows inh reverse xS diffClass sameCell r=.058 vs diffCell .083 (no cell effect cross-study);
wS diffClass sameCell .41 vs diffCell .31 (within-study; shared-control confound possible). Declarations below are made with that knowledge.

## Units
Experiment = study x cell x agent (audit_I19 inh set, binder-masked) or study x cell x stimulus class (stim set). Reverse = up for inhibitors, down for stimuli.
Profile = per S/T Cantley kinase, enrichment of motif (pct>=90) among reverse responders vs non-responders (I05 definition) and a continuous variant
(Spearman of motif pct vs reverse magnitude). Profiles centred per experiment (shape), then per kinase across experiments (remove the global reverse signature).

## H1 (cell-intrinsic repertoire)
Statistic: mean profile r for pairs (same cell, different class) minus (different cell, different class), with the pair type restricted
(a) cross-study (independent controls; primary test), (b) within-study multi-cell studies (Hijazi 3 cells, decryptM), with shared-control correction.
Null: block permutation of cell labels across study x cell units (1000 perms); CI: study cluster bootstrap.
Declared prediction if H1 true: delta>0, perm p<0.05, and delta(reverse) > delta(primary). Declared falsification: cross-study delta <= 0 or CI includes 0
while within-study delta disappears after removing the shared control component -> verdict "not supported / artefact".

## H2 (baseline expression predicts repertoire)
Statistic: per cell line, Spearman over kinases between cell-relative expression (z across cell lines, per kinase) and cell-relative reverse enrichment
(z across cell lines, per kinase); mean over cells vs permuted cell labels. Within-study version (Hijazi MCF7/HL60/NTERA2; decryptM cells) is primary.
Pre-declared kinase set: S/T Cantley kinases with expression data; secondary: family-aggregated (max expression in motif family/cluster).
Declared prediction: mean matched r > 0, perm p < 0.05; also specific hypothesis: MAP4K4/TNIK/MINK1/NEK expression vs their reverse LOR, positive.
Falsification: matched r <= permuted.

## H3 (alternatives)
Coverage (sites measured; restrict profiles to sites measured in both experiments of a pair), study (cross-study tests; studies with >=2 cells),
baseline phosphorylation/saturation (no absolute baseline in compendium -> proxy: site's typical level not available; use protein abundance / site detection frequency).
Declared: if cell effect (H1) survives common-site restriction and cross-study comparison, coverage/study do not explain it.

## H4 (optional)
Chronos dependency on kinase k in cell c (z across cells) vs reverse enrichment of k in c: same matched-vs-permuted statistic. Prediction: negative relation (more dependent -> more compensatory reverse); weak prior.

## Results log (after declaration)
- Downloads: DepMap 24Q2 (figshare 25880521: Model.csv, OmicsExpressionProteinCodingGenesTPMLogp1) + Gygi CCLE proteome OK in ~1 min (DepMap portal API blocked by Cloudflare).
  47/73 compendium cell labels mapped to DepMap RNA, 26 to proteome, 38 to Chronos 26Q1; HEK293/RPE1 only in HPA (8 cells); WIDR->HT29 (STR identical).
- H1 motif level (T1; columns r_same, r_diff, delta — an earlier quick read mixed columns, corrected here): cross-study diff-class delta reverse -0.03..+0.12 (main A dpct raw -0.012 [-0.20,0.19] p=.62; inh-only dpct +0.09/+0.12 p=.12/.055 = weak trend; lor90 -0.02/+0.03), primary -0.10..-0.03. NOT SUPPORTED.
  Within Hijazi (shared DMSO per cell): reverse clusters by cell (same cell diff class r .15-.18) more than by drug (same drug diff cell .06-.09); primary by drug (.15-.21 vs .11-.17).
  Survives per-cell median centring and restriction to 4,908 sites measured in all three cells (T4) -> not shared-control additive artefact, not coverage.
  decryptM (independent per-drug controls; only 11-48 same-cell pairs): reverse same-cell excess -0.11..+0.08, primary +0.18..+0.30 (consistently larger).
  Within-study permutation is uninformative (same cell == same unit).
- H1 site level (T3/T4): cross-study diff class inh reverse delta LOR +0.20 [-0.29,0.54] p=.20; same class +0.03 (primary +0.30). MCF7 Hijazi vs MCF7 Ochoa (2106 pairs): reverse LOR -0.10 same cell vs -0.15 diff -> reverse sites do not replicate even in the same cell line at the same time point.
- Interpretation: "cell-line specificity" of the reverse layer = unit (study x cell culture) specificity, not genotype-intrinsic.
- H2 (T5, T6, T8): all matched-vs-permuted tests null (RNA, proteome, Chronos; per kinase or motif family; all units or within study): |mean rho|<=0.035, p>=0.16.
  Targeted 132 tests: 4 at p<.05 (< 6.6 expected), signs inconsistent (MINK1 RNA rho -0.31). Within Hijazi: reproducible cell-specific reverse profiles (split-half .41-.94) do not track own RNA (diag rho -0.07..+0.04).
  Positive control (T7): phosphosite detection per gene vs RNA, matched rho .135 (all genes)/.162 (kinases) vs null 0 +- .01, p=.003 -> mapping & test have power.
- H4: null (Chronos rows of T5; MAP4K4 rho +0.42 p=.01 in opposite direction to prediction, isolated).
