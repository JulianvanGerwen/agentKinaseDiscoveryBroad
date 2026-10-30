# I26: auditing PSP kinase–substrate edges with perturbation data (primary-layer decreases only)
**Setup.** 29 inhibitor-resolvable kinase groups (family level, e.g. ERK1/2, AKT1-3, CDK1/2/5/9) built from I03 LoF conditions (461 conditions). Edges: 9,305 PSP edges in the infra psp_z library, 5,688 with the site in the compendium. Evidence is the per-study mean robust z over LoF conditions; decreases only. **Positive control, within the same experiment:** a condition counts as "worked" when the group's established PSP substrates fall (in-vivo edges, site LT_LIT≥3; Stouffer ≤−2.5, ≥3 sites, leaving out the edge under test). 146/476 conditions pass, from 23 studies. Classes need ≥2 studies with verified inhibition: **supported** = ≥2 studies and ≥50% with z≤−1; **unsupported** = <25% of studies with z≤−1 and mean >−0.5; otherwise **ambiguous**.
**1. Coverage and rates.** Testable: 680 kinase-edges / 544 (group, site) units in 15 groups; 93% of edges cannot be tested. Most untestable sites were measured in <2 LoF studies (8,017); for 608 the inhibition was never verified twice. Groups with no verified inhibition: CK2 (CX-4945 does not move canonical CK2 sites at 1 h), p38, JAK, PLK1, TAK1, CDK7/12/13, FLT3. Units: **supported 17% (93), unsupported 50% (270), ambiguous 33%**. By group, supported vs matched random sites: CDK1/2/5/9 12% vs 6%, ERK 16% vs 4.5%, mTOR 26% vs 6%, AKT 17% vs 1.5%, EGFR/ERBB2 53% vs 7%.
**Calibration (task 4).** Random relabelled control: same group, random measured sites matched on residue class and number of studies, ×20. Supported **5.3%**, unsupported **66%**. So PSP edges are enriched about 3× for support, but only ~12% (supported) and ~25% (excess non-unsupported) of testable edges carry signal beyond random. **Gold edges** (LOO) are themselves only 22% supported and 44% unsupported. The assay is therefore insensitive, and "unsupported" ≠ "wrong". The ordering PSP > random holds across all 24 threshold settings (control Stouffer none/−1.5/−2.5/−4 × decrease −0.5/−1/−1.5 × gold LT≥1/3; fig1C). The positive-control filter changes rates little (supported 14% with no control vs 17%).
**Relation to PSP annotation** (logit, group fixed effects, n=518):
- More low-throughput references: support OR 1.23/log-unit [1.02, 1.48]; mean-z β −0.19 [−0.30, −0.08].
- Cantley motif percentile: weak (β −0.009 z/pct, p=0.02; no effect on the supported class).
- In-vivo vs in-vitro-only flag: no difference (OR 0.87 / 1.51, CIs span 1). Supported is 16% for vitro-only vs 18% for vivo+vitro.
- Kinase family: CMGC (CDK/ERK) has the lowest support (13%), TK/STE the highest. Those groups are few.
**2. Corrections** (`tables/corrections*.tsv`). 91/270 unsupported units have an alternative group that is supported, 54 of them with verified inhibition. Filtered for verified support, ≤2 groups supporting the site, alt z≤−1.5 and motif ≥85 (alt group or its canonical downstream kinase): **27 edges**. Of these, 18 are novel alternatives, 7 have the alternative already in PSP, and 2 are pathway-consistent. Dominant pattern: AKT/CDK/PKC-annotated basophilic or SP sites that behave as **MEK–ERK–RSK** sites.
Literature check (PSP KS file has no PMIDs; checked via web):
- **EPHA2 S897** (AKT→RSK): confirmed by MEK/RSK-but-not-AKT inhibitor studies.
- **CHEK1 S280** (AKT→RSK): Li 2012 shows RSK, not AKT.
- **CTTN S405** (PKCδ/CDK→ERK): classic ERK site.
- **MAP2K2 T394** (CDK5→ERK feedback).
- **DNMT1 S143**: AKT evidence was in-vitro kinase assay plus knockdown.
- **PIK3C2A S259**: CDK1 inferred from G2/M correlation.
- **FOXK2 S373**: the CDK paper mapped S368/S423; our data show mTOR dependence.
- **EDC3 S161**: AKT vs PIM conflict; our data show neither AKT nor mTOR, but RSK-pathway dependence.
- Context, not error: GRK2 S670 (CDK2 G2/M vs ERK feedback), CUEDC2 S110 (mitotic CDK1).
- No literature found for DBNL S269, TFPT S180, AKAP11 T1100. These are new RSK/ERK candidates.
**3. Activity inference** (PSP Stouffer, ≥3 substrates; 394 inhibitor conditions, 32 studies). The audit is recomputed without the evaluated study (LOSO) and in 20 random half-splits. Change in AUROC (Δ):

| Variant | LOSO, all conditions | LOSO, conditions with an edit | Half-splits, conditions with an edit |
|---|---|---|---|
| Drop unsupported | **+0.005** (studies p=0.005) | +0.012 | +0.012 |
| Random drop, same number, non-supported pool | +0.003 | +0.008 | +0.009 |
| Random drop, any edge | −0.002 | −0.004 | −0.002 |
| Drop unsupported + ambiguous | +0.009 | | |
| Adding corrections | −0.001 | | |

The gain is real but small, and mostly comes from removing non-responsive edges, whichever ones are chosen. core-v3 was not re-run.
**Novelty.** Hijazi 2020 used inhibitor data to *score* activity, not to audit edges. decryptM 2.0 (Bayer 2025 preprint; ~31% of PSP pairs potency-coherent) audits by dose coherence within one lab. benchmarKIN ranks libraries. New here: a cross-lab audit with within-experiment positive controls; a random-site calibration showing even gold edges are mostly "unsupported" (an upper bound on what any acute-perturbation audit can claim); and a ranked, literature-checked correction list. Our ~17% supported / ~33–50% non-unsupported range brackets decryptM's 31%.
**Negatives and caveats.**
- 93% of edges are untestable.
- Results are dominated by Hijazi and decryptM: excluding Hijazi leaves 61 unsupported units.
- Alternatives are mostly pathway-level (MEK/RAF/ALK inhibitors). The direct kinase is inferred by motif.
- Acute (≤1–2 h), asynchronous cells: mitotic CDK1 edges and slow-turnover CK2 sites are systematically unassessable, which gives a context bias.
- The gold set is PSP-derived (LOO).
- Upstream and downstream confounds remain.
**Self-score:** Compelling 3, Novel 3, Promising 3. The correction list (EPHA2, CHEK1, DBNL, MAP2K2) is the most actionable output. The headline caveat is the insensitivity shown by the gold edges.
Files: `tables/` (edge_audit, per_group/kinase_support, random_control, sensitivity_grid, relations_*, corrections_curated, activity_benchmark_*), `figures/fig1–3`, `scripts/b_I26/run_all.sh`.
