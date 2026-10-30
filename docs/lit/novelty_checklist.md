# Novelty checklist: known findings about kinase signaling from phosphoproteomics

Use this table to judge whether one of our results is novel.

**Status values:**
- **Known**: a published claim at the same level of specificity.
- **Partly**: a related claim exists, but our angle (scale, context, quantity) may still be new.
- **Open**: no prior claim found. Re-check before publishing.

**Evidence tags:**
- **[V]**: claim checked against the abstract or full text this session.
- **[M]**: from memory.

Full context is in `/workspace/LITERATURE.md`.

## A. Annotation, inference and benchmarks
| # | Finding | Status | Ref |
|---|---|---|---|
| A1 | <5% of phosphosites have a known kinase; ~90% of annotated sites belong to the top 20% of kinases; >150 kinases have no known substrate | Known [V] | Needham 2019, PMID 30670635 |
| A2 | The number of known substrates per kinase does not correlate with disease relevance | Known [V] | Needham 2019 |
| A3 | Z-test ≈ GSEA > KS-test/Wilcoxon/MLR for activity inference; in vivo > in vitro > in silico substrates; ≤5 substrates gives poor inference | Known [V] | Hernandez-Armenta 2017, PMID 28200105 |
| A4 | The KS library matters more than the algorithm; combined curated libraries are best; NetworKIN helps on tumours; Kinase Library top-15 targets hurt | Known [V] | Müller-Dott 2025, PMID 40404650 |
| A5 | Network propagation (RoKAI) improves activity inference and robustness to missing substrates | Known [V] | Yılmaz 2021, PMID 33608514 |
| A6 | Motif-only (Kinase Library) activity inference is competitive with curated-library methods | Known [V] | Lussana 2024 (PhosX), PMID 39563468 |
| A7 | Study-bias-matched nulls are needed for patient-specific activity | Known [V] | Crowl 2022 (KSTAR), PMID 35879309 |
| A8 | Literature networks beat data-derived KS networks for EGF pathway recovery; KS networks cover <5% of sites | Known [V] | Garrido-Rodriguez 2026, PMID 41741430 |
| A9 | Intrinsic substrate specificity atlases exist for 303 S/T kinases and the full Tyr kinome; many S/T kinases prefer primed (phospho) residues | Known [V] | Johnson 2023, PMID 36631611; Yaron-Barir 2024, PMID 38720073 |
| A10 | In vitro lysate screens give about 175k KS pairs for 385 kinases, with low overlap with PSP | Known [V] | Sugiyama 2019, PMID 31324866 |
| A11 | ATP-competitive inhibitors can increase phosphorylation of their target's own activation loop or hydrophobic motif (AKT) | Known [V] | Okuzumi 2009, PMID 19465931 |
| A12 | A systematic, kinome-wide catalogue of inhibitor-induced paradoxical target hyperphosphorylation | Open [M] | none found |
| A13 | Normalising sites to protein abundance lowers tumour-benchmark performance | Partly [V summary] | bioRxiv 2024.06.27.601117 |
| A14 | Phosphatase activity inference from phosphoproteomics (benchmarked) | Open [M] | none found |
| A15 | Data-driven, cross-study consensus "perturbation footprints" (PROGENy-style) for kinases from phosphoproteomics | Open/Partly [M] | PTMsigDB PERT sets (Krug 2019, PMID 30563849) are per-study signatures |

## B. Compendia, function, evolution, networks
| # | Finding | Status | Ref |
|---|---|---|---|
| B1 | Activity atlas for 215 kinases × 399 conditions; frequently regulated kinases are network-central; co-regulation predicts complexes and KS; activity profiles act as state fingerprints | Known [V] | Ochoa 2016, PMID 27909043 |
| B2 | New AKT1 targets, stem-cell differentiation mediators and Salmonella modulators from the kinase atlas | Known [V] | Ochoa 2016 |
| B3 | Signed kinase→kinase regulatory network predicted kinome-wide; denser than textbook pathways | Known [V] | Invergo 2020, PMID 32437683 |
| B4 | Inhibitor phosphoproteomics (61 drugs) infers >1,500 kinase–kinase relations; PIK3CA-WT breast cells use MAPK-dependent circuitry; TTK important in AML; kinase expression is only partly coupled to activity | Known [V] | Hijazi 2020, PMID 31959955 |
| B5 | Community detection in empirical kinase networks proposes new pathway members | Known [V title] | Basanta 2023, PMID 37352361 |
| B6 | Functional score for 119,809 human sites from 59 features; SMARCC2 sites in neuronal differentiation | Known [V] | Ochoa 2020, PMID 31819260 |
| B7 | Most phosphosites are weakly constrained, and many may be non-functional | Known [V] | Landry 2009, PMID 19349092; Lienhard 2008, PMID 18603430 |
| B8 | PTM hotspots in domain families mark functional regions; only a fraction of sites are functional | Known [V] | Beltrao 2012, PMID 22817900; Strumillo 2019, PMID 31036831 |
| B9 | Phosphosites turn over rapidly in evolution; ancient sites are enriched at interfaces and more functional | Known [V] | Studer 2016, PMID 27738172 |
| B10 | Kinase motifs largely date to the last eukaryotic common ancestor; recent duplicates keep specificity | Known [V] | Bradley 2019, PMID 31233486 |
| B11 | CK2 site "substrate quality" determines function and evolution | Known [V title] | Bradley 2024, PMID 38861992 |
| B12 | Positive selection for tyrosine loss in metazoans | Known [V title] | Tan 2009, PMID 19589966 |
| B13 | Tyr and S/T phosphorylation have distinct regulatory natures (deep HeLa phosphoproteome) | Known [V title] | Sharma 2014, PMID 25159151 |
| B14 | Widespread full site occupancy in mitosis | Known [V title] | Olsen 2010, PMID 20068231 |
| B15 | Kinase AP-MS interactome of >300 kinases with >5,000 new interactions | Known [V] | Buljan 2020, PMID 32707033 |
| B16 | Node abundance shapes signaling network state (mass cytometry) | Known [V title] | Lun 2017, PMID 28092656 |
| B17 | Site age or conservation vs *dynamic regulatory behaviour* (breadth, dose sensitivity) in human perturbation compendia | Open [M] | none found |
| B18 | ICA or decomposition-based site "phosphomodules" across a large perturbation compendium, transferred to tumours | Open [M] | Ochoa 2016 is kinase-level; CoPheeMap is a tumour co-regulation network |
| B19 | Dark-kinase activity atlas across thousands of conditions | Open/Partly [M] | Ochoa 2016 covered 215 mostly well-studied kinases |
| B20 | Mutations at or near phosphosites rewire signaling (ActiveDriver, MIMP, network-attacking mutations, pY switches) | Known [V] | Reimand 2013, PMID 24089029; Wagih 2015, PMID 25938373; Creixell 2015, PMID 26388441; Lundby 2019, PMID 31585087 |
| B21 | Genetic variants affecting phosphorylation differ from those affecting abundance; phospho traits correlate better with physiology (yeast phosphoQTL) | Known [V] | Grossbach 2022, PMID 35574625 |

## C. CPTAC tumour phosphoproteomics
| # | Finding | Status | Ref |
|---|---|---|---|
| C1 | BRCA: amplicon-associated hyperphosphorylated kinases CDK12, PAK1, PTK2, RIPK2, TLK2 (besides ERBB2) | Known [V] | Mertins 2016, PMID 27251275 |
| C2 | BRCA: 5q loss (CETN3, SKP1) → raised EGFR; SKP1 loss → SRC | Known [V] | Mertins 2016 |
| C3 | BRCA: CDK4/6 activity (PTM-SEA) tracks proliferation in HR+ but not TNBC; Rb phosphorylation for CDK4/6i; mTOR activity higher in Rb-high TNBC | Known [V] | Krug 2020, PMID 33212010 |
| C4 | COAD: Rb phosphorylation is linked to proliferation and reduced apoptosis (CDK targeting rationale) | Known [V] | Vasaikar 2019, PMID 31031003 |
| C5 | ccRCC: CDK1 and MAPK1 substrates are top tumour-vs-NAT events; WEE1→CDK1 and CDK7→MCM2 are up; CDK1 pY15 is higher in aggressive tumours; EGFR is broadly up | Known [V] | Clark 2019, PMID 31675502 |
| C6 | ccRCC: EIF4EBP1 is a MAPK/mTOR convergence node, suggesting combined mTOR+MEK | Known [V] | Clark 2019 |
| C7 | ccRCC: STAT1 phosphorylation in the CD8-inflamed immune subtype | Known [V] | Clark 2019 |
| C8 | UCEC: CDK4 activity is higher in CNV-high endometrioid tumours; CDK12 is hyperphosphorylated | Known [V] | Dou 2020, PMID 32059776 |
| C9 | UCEC: PIK3R1 in-frame indels → raised AKT phosphorylation and AKTi sensitivity | Known [V] | Dou 2023, PMID 37567170 |
| C10 | UCEC: CTNNB1 hotspot mutations cluster at the S45 degron phosphosites | Known [V] | Dou 2023 |
| C11 | LUAD: ALK fusions → outlier ALK pY1507 | Known [V] | Gillette 2020, PMID 32649874 |
| C12 | LUAD: EGFR-mutant → PTPN11 pY62 (protein unchanged); ALK-fusion → PTPN11 pY546/pY584 | Known [V] | Gillette 2020 |
| C13 | LUAD: KRAS-mutant → SOS1 pS1161, MAPK1 and TAK1 activity; EGFR-mutant → ROCK1 | Known [V] | Gillette 2020 |
| C14 | LUAD: kinase phospho-outliers PRKCD (KRAS-mut), BRAF (TP53-mut), WEE1 (EML4-ALK), matched to shRNA/CRISPR dependency | Known [V] | Gillette 2020 |
| C15 | LUAD: STK11 mutation → immune-cold; STK11 → SMAD4 pS138 | Known [V] | Gillette 2020 |
| C16 | GBM: PTPN11 and PLCG1 phosphorylation as RTK-pathway switches; PTPN11 central across high-grade glioma | Known [V] | Wang 2021, PMID 33577785; Liu 2024, PMID 38981438 |
| C17 | HNSCC: two modes of EGFR activation; ligand abundance predicts anti-EGFR mAb benefit | Known [V] | Huang 2021, PMID 33417831 |
| C18 | HNSCC: CIN subtype has high CDK4/6 activity (Rb hyperphosphorylation); CHEK2 activating site hyperphosphorylated | Known [V] | Huang 2021 |
| C19 | LSCC: CDK5 activity is up in the Classical/NRF2 subtype; NSD3 in FGFR1-amp tumours; RB1 phosphorylation informs CDK4/6i | Known [V] | Satpathy 2021, PMID 34358469 |
| C20 | OV: mitotic kinases and replication stress implicated in CIN | Known [V] | McDermott 2020, PMID 32529193 |
| C21 | Pan-cancer: HRD/hypoxia tumours have raised ATM, CHEK1, CHEK2 and CDK1/2/4/6 activity (PTM-SEA); PARP1 S179/S782 and EXO1 S714 up | Known [V] | Geffen 2023, PMID 37582358 |
| C22 | Pan-cancer: MMRd tumours with ATM indels show reduced ATM→PRKDC S3205 and PNKP S114/T118 | Known [V] | Geffen 2023 |
| C23 | Proximal lysine acetylation (+3K) blocks CDK1 phosphorylation (RSF1 S1375); acetylation alters S/T kinase specificity | Known [V] | Geffen 2023 |
| C24 | Pan-cancer: most driver genes converge on similar sequence-based kinase-activity states | Known [V] | Li 2023 Cell, PMID 37582357 |
| C25 | Pan-cancer immune subtypes: FAK (PTK2), CDK19 and CDK20 more active in hot tumour cells; FYN/LYN/LCK lower in cold tumours; CDK1/2 high in CD8− subtypes | Known [V] | Petralia 2024, PMID 38359819 |
| C26 | ATR kinase activity is driven by its own CNV in HNSCC | Known [V] | Petralia 2024 |
| C27 | Protein overexpression or hyperactivation in CPTAC + DepMap → druggable dependencies; proteogenomic synthetic lethality with TSG loss | Known [V] | Savage 2024, PMID 38917788 |
| C28 | CPTAC co-phosphorylation network (26k sites) predicts 24k KS links | Known [V] | Jiang 2025 (CoPheeKSA), Nat Commun 16:2766 |
| C29 | Kinase protein extremes and activating sites as tumour truth for activity benchmarks; COAD/OV/PDAC excluded for low site–protein correlation | Known [V] | Müller-Dott 2025 |
| C30 | Driver-gene cis-effects (mutation/CNV raising own protein and phosphosite) are pervasive (e.g., TP53, EGFR up; STK11, KEAP1, NF1, RB1 down) | Known [V] | Gillette 2020; Li 2023 |
| C31 | Kinase-pair or network-state (not single-kinase) differences between tumours and perturbation-derived wiring | Open [M] | none found |
| C32 | MR-style (cis-CNV instrument) causal kinase→substrate effect estimates kinome-wide | Open/Partly [M] | Petralia 2024 (descriptive, one kinase) |
| C33 | Phospho-CMap: scoring tumours with proteome-scale inhibitor footprints | Open/Partly [M] | P100 CMap uses ~96 sites in cell lines (Abelin 2016, PMID 26912667) |
| C34 | Pareto/archetype geometry of tumour kinase activity | Open [M] | none found |

## D. Inhibitors, feedback, resistance, off-targets
| # | Finding | Status | Ref |
|---|---|---|---|
| D1 | MEK inhibition → ERK loss → c-MYC degradation → RTK induction (TNBC) | Known [V] | Duncan 2012, PMID 22500798 |
| D2 | AKT inhibition relieves feedback suppression of RTK expression and activity | Known [V] | Chandarlapaty 2011, PMID 21215704 |
| D3 | mTOR kinase inhibition gives biphasic, feedback-dependent AKT regulation | Known [V] | Rodrik-Outmezguine 2011, PMID 22140653 |
| D4 | RAF inhibitors relieve ERK feedback in BRAF-V600E cells; paradoxically activate RAF dimers in BRAF-WT | Known [V] | Lito 2012, PMID 23153539; Poulikakos 2010, PMID 20179705; Hatzivassiliou 2010, PMID 20130576 |
| D5 | PI3Kβ inhibition relieves feedback on PI3Kα in PTEN-mutant tumours | Known [V] | Schwartz 2015, PMID 25544636 |
| D6 | Lapatinib kinome reprogramming is blocked by BET inhibitors; lapatinib resistance involves glycolysis phospho-reprogramming | Known [V] | Stuhlmiller 2015, PMID 25865888; Ruprecht 2017, PMID 28209619 |
| D7 | Cells resistant to the same PI3K/mTOR inhibitor evolve heterogeneous network states that predict combination response | Known [V] | Wilkes 2015, PMID 26060313 |
| D8 | BRAF-V600E CRC intrinsic resistance via PDPK1 and PRKCA; RPS6KB1/PIM1 vulnerabilities in melanoma | Known [V] | Coppé 2019, PMID 31160710 |
| D9 | Dose- and time-resolved PTM profiling: 31 drugs × 13 lines, 1.8M curves; rituximab kills by BCR over-activation | Known [V] | Zecha 2023, PMID 36926954 |
| D10 | Kinobeads target landscapes: 243 clinical drugs; 1,183 compounds with selective probes for 72 kinases | Known [V] | Klaeger 2017, PMID 29191878; Reinecke 2024, PMID 37904048 |
| D11 | Phosphoproteomic time series + activity inference reveal adaptive rewiring to clinical inhibitors in CRC (VESPA) | Known [V] | Rosenberger 2023, PMID 36824919 |
| D12 | CDKs play a key role in TNF signaling | Known [V title] | Tanzer 2021, PMID 34663829 |
| D13 | Causal signaling edges are context (cell-line/stimulus)-specific (RPPA) | Known [V title] | Hill 2017, PMID 28017544 |
| D14 | Proteome-wide EC50 vs Kd_app comparison to quantify amplification or ultrasensitivity and network depth | Open [M] | Theory: Huang & Ferrell 1996, PMID 8816754 |
| D15 | Functional (cellular) off-target calls from site dose–response matching non-target Kd_app | Open/Partly [M] | Engagement examples in Zecha 2023 and Reinecke 2024 |
| D16 | Cross-cell-line invariance atlas of feedback (sites rising upon inhibition) | Open [M] | Axis-specific studies only |
| D17 | Kinase activity scores from cell-line phosphoproteomes correlate with drug sensitivity | Known [V] | Koenig 2026, PMID 42816593; Frejno 2020, PMID 32686665 |

## E. Dynamics and network inference
| # | Finding | Status | Ref |
|---|---|---|---|
| E1 | Global site-specific EGF phosphorylation dynamics | Known [V title] | Olsen 2006, PMID 17081983 |
| E2 | Seconds-scale EGFR signaling dynamics | Known [V title] | Reddy 2016, PMID 26929352 |
| E3 | AKT directly regulates mTORC2 (insulin time course) | Known [V title] | Humphrey 2013, PMID 23684622 |
| E4 | Early responses are specific, later ones promiscuous (yeast) | Known [V title] | Kanshin 2015, PMID 25704821 |
| E5 | Data-only causal network inference from phospho time courses is weak; prior knowledge helps (HPN-DREAM) | Known [V] | Hill 2016, PMID 26901648 |
| E6 | Prior-knowledge ILP/causal methods (PHONEMeS, CARNIVAL, COSMOS, CausalPath) reconstruct pathways from phospho data | Known [V] | PMIDs 26354681, 31728204, 33502086, 34179843 |
| E7 | Inter-individual human phosphoproteomes (insulin resistance) reveal selective mTORC1 defects; MINDY1 S441 is an AKT substrate | Known [V] | Needham 2024, PMID 39577414 |
| E8 | Cross-study human meta-analysis of temporal response order (early robust sites enriched for direct KS motifs) | Open [M] | none found |
| E9 | Modular response analysis or local response coefficients applied to proteome-scale multi-inhibitor data | Open [M] | Theory: Kholodenko 2002, PMID 12242336 |
| E10 | Nested-effects or footprint-containment hierarchy of kinases from inhibitor panels | Open [M] | Theory: Markowetz 2005, PMID 16159925 |
