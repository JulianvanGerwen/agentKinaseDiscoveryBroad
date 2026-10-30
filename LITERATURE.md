# LITERATURE.md: What phosphoproteomics has taught us about kinase signaling, and what is still open

Compiled 2026-10-29 for the kinase-signaling discovery project (see PROMPT.md, RESEARCH_TREE.md).
Companion files:
- `docs/lit/novelty_checklist.md`: a table of known findings, used to judge whether our results are novel.
- `docs/lit/SUMMARY.md`: a one-page digest.

**Scope.** This review covers what has been *learnt about kinase signaling* from phosphoproteomics data. It does not re-review activity-inference methods in depth. Those are covered in two read-only files:
- `kinaseActivityMain/LITERATURE.md`: methods, KS libraries and benchmarks.
- `kinaseSubMain/LITERATURE.md`: substrate prediction.

Section 1 only summarises their conclusions.

**Evidence tags:**
- **[V]**: the PubMed record, abstract or full-text sentence was checked this session. PMIDs/DOIs were checked through NCBI E-utilities.
- **[M]**: from memory or secondary sources and not re-checked. Treat specifics as unverified.

A PMID with no tag means the PMID and title were verified, but the specific claim may be from memory where marked.

---

## 0. Ten-line orientation

1. **Coverage, not the algorithm, limits inference.**
   - Fewer than 5% of phosphosites have a known kinase. About 90% of annotated sites belong to the top 20% of kinases, and more than 150 kinases have no known substrate (Needham 2019, PMID 30670635 [V]).
   - Benchmarks agree that the library matters more than the statistic (Hernandez-Armenta 2017, PMID 28200105; Müller-Dott 2025, PMID 40404650 [V]).
2. **Compendia exist and have been mined, at kinase level, for these questions:**
   - kinase co-regulation (Ochoa 2016, 215 kinases × 399 conditions, PMID 27909043 [V])
   - signed kinase–kinase circuits (Invergo 2020, PMID 32437683 [V]; Hijazi 2020, PMID 31959955 [V])
   - site function (Ochoa 2020, 119,809 sites, PMID 31819260 [V])

   What has *not* been mined at scale is site-level modular structure, dose–response structure, or cross-context invariance.
3. **CPTAC has made many specific kinase claims.** Most are driver→phosphosite cis/trans effects, "hyperphosphorylated kinase outliers", CDK4/6/Rb, DNA-damage kinases, EGFR/PTPN11 and immune-subtype kinases (§3; checklist). Generic claims such as "CDK activity is high in proliferative tumours" or "EGFR-mutant tumours have high EGFR/PTPN11 phosphorylation" are *not* novel.
4. **Feedback and rewiring after inhibition are known for a handful of axes**: MEK→MYC→RTK, AKT→FOXO→RTK, RAF paradox, mTORC1→IRS1/AKT. No proteome-scale, cross-context catalogue of feedback responses with dose resolution exists. decryptM and Hijazi make one possible (§4).
5. **Causal and temporal inference** relies on prior-knowledge-constrained optimisation (PHONEMeS, CARNIVAL, COSMOS, CausalPath). The DREAM challenge showed data-only network inference is weak (Hill 2016, PMID 26901648 [V]). Literature networks still win for EGF pathway recovery (Garrido-Rodriguez 2026, PMID 41741430 [V]).

---

## 1. Kinase activity inference: methods and what large-scale applications found

### 1.1 Methods (condensed; details in kinaseActivityMain/LITERATURE.md)

**Footprint/enrichment methods** use known kinase–substrate (KS) sets:
- KSEA, a z-test (Casado 2013, PMID 23532336 [V]; KSEAapp: Wiredja 2017, PMID 28655153 [V]).
- GSEA-type methods.
- PTM-SEA/PTMsigDB, with directional site signatures (Krug 2019, PMID 30563849 [V]).
- KinSwingR, a PWM-network "swing" score (Engholm-Keller 2019, PLoS Biol [M]).
- KARP (Wilkes 2017, PMID 28674151 [V]).
- INKA, which combines A-loop phosphorylation with substrate evidence (Beekhof 2019, MSB [M]).
- KEA3 (gene-level ORA).
- decoupler ulm/mlm/viper (Badia-i-Mompel 2022, PMID 36699385 [V]).

**Network and latent-variable methods:**
- RoKAI: network-propagated site values (Yılmaz 2021, PMID 33608514 [V]).
- KSTAR: binarised NetworKIN graphs with study-bias-matched nulls (Crowl 2022, PMID 35879309 [V]; v1.2 2026 preprint).
- IKAP: joint activity × affinity fit (Mischnik 2016, PMID 26628587 [V]).
- JUMPsem: per-kinase confirmatory factor analysis (Kong 2025, PMID 39838083 [V]).
- PhosR: kinase–substrate scoring from temporal co-regulation (Kim 2021, PMID 33626354 [V]).

**Motif-only ("annotation-free") methods:**
- **Kinase Library:** positional-scanning specificity for 303 S/T kinases (Johnson 2023, PMID 36631611 [V]) and the tyrosine kinome (Yaron-Barir 2024, PMID 38720073 [V]). Enrichment uses each site's top-15 predicted kinases.
- **PhosX:** a Kinase-Library-based GSEA (Lussana 2024, PMID 39563468 [V]).
- **Kinex** (bioRxiv 2023 [V in the activity review]).

**Benchmarks:**
- Hernandez-Armenta 2017 (PMID 28200105 [V]) found:
  - Z-test ≈ GSEA, both better than KS-test, Wilcoxon and ridge MLR. Median AUC was about 0.72–0.73.
  - in-vivo substrates > in-vitro > NetworKIN.
  - Kinases with ≤5 quantified substrates were poorly scored.
- benchmarKIN (Müller-Dott 2025, PMID 40404650 [V]) used 19 methods × 8 libraries on perturbation and CPTAC tumour benchmarks:
  - Most methods were similar. The z-score was near-best.
  - Combined curated libraries were best.
  - NetworKIN helped on tumours.
  - Naive Kinase Library or Phosformer top-15 targets *hurt* (from the activity review [V]).
- Garrido-Rodriguez 2026 (PMID 41741430 [V]) built an EGF benchmark with new deep time courses. The literature network recovered EGF pathway interactions best, and KS networks covered <5% of sites [V in activity review].

### 1.2 What large-scale applications concluded about biology (not methods)

- **Kinase regulation is low-dimensional and hub-dominated.**
  - Across 399 conditions, the most frequently regulated kinases are central in the signaling network, and co-regulated kinase pairs predict complexes and KS relations (Ochoa 2016, PMID 27909043 [V]).
  - Kinase-activity profiles act as fingerprints that match related and opposing cell states. They identified mediators of stem-cell differentiation, modulators of Salmonella infection and new AKT1 targets [V abstract].
- **Inhibitor phosphoproteomics can map kinase–kinase topology.**
  - Hijazi 2020 (PMID 31959955 [V]) used 61 inhibitors × 3 cell lines to derive more than 6,000 "putative downstream target" sites and more than 1,500 kinase–kinase relations.
  - Applied to AML and breast tumours, it showed that kinase expression and activity are only partially coupled.
  - It showed that PIK3CA-wild-type breast cells adopt MAPK-dependent circuitry, and that TTK matters in AML [V].
  - Follow-up: community detection in these empirical networks proposes new pathway members (Basanta 2023, PMID 37352361 [V title]).
- **Activity from single samples is feasible but biased toward well-studied kinases.**
  - KARP, INKA and KSTAR were developed for patient samples.
  - KSTAR emphasised that "patient-specific kinase activities" need study-bias-matched nulls (PMID 35879309 [V]).
- **Motif-only inference is competitive for some kinases** (PhosX [V]). But family members are indistinguishable at the motif level. This is a hard ceiling for kinase-specific (vs family-specific) claims.
- **Deep cell-line atlases:**
  - Frejno 2020 (PMID 32686665 [V]): 125 cancer cell lines, 55k sites. It related phosphosites to drug response; examples are PGR phosphorylation and raloxifene sensitivity, and AK1 and cytarabine.
  - Koenig C et al. 2026 (NSMB, PMID 42816593 [V]; authors corrected in the briefing fact-check): 33 cell lines, more than 200k class-I sites. A combined proteome+phosphoproteome kinase activity score ("Cscore") gave cell-line-specific kinase vulnerabilities that correlate with drug sensitivity. **This is a direct novelty threat for "kinase activity predicts drug sensitivity in cell lines" claims.**
- **Human in-vivo, inter-individual variation:**
  - Needham 2024 (Cell Metab, PMID 39577414 [V]): 114 muscle biopsies from 19 men, with insulin and exercise.
  - Insulin resistance selectively reduces mTORC1 and non-canonical insulin signaling, and exercise "primes" part of the network.
  - It identified MINDY1 S441 as an AKT substrate.
  - This shows *between-person* variation in phosphoproteomes can reveal kinase biology, a natural-experiment design.

### 1.3 Known limitations that create opportunities

1. **Study bias.** The √m scaling of z-scores favours well-annotated kinases. Benchmarks are dominated by about 25–30 kinases (RoKAI: median 29 vs 4 substrates [V in activity review]).
2. **Paralog indistinguishability.** This applies to motif-based inference and to inhibitor truth sets.
3. **The activation-loop paradox.** ATP-competitive inhibitors can *increase* phosphorylation of their target's activation loop or hydrophobic motif (Okuzumi 2009, PMID 19465931 [V]). RAF inhibitors paradoxically activate RAF dimers (Poulikakos 2010, PMID 20179705 [V]; Hatzivassiliou 2010, PMID 20130576 [V]).
4. **Protein abundance confounds tumour phosphosite variation.** In the tumour benchmark, normalising sites to protein *lowered* performance (bioRxiv 2024.06.27.601117 [V via search summary]). This suggests some "activity" in tumours is abundance-driven.
5. **Phosphatases are almost absent** from activity inference. No benchmarked phosphatase footprint resource was found.

---

## 2. Large compendium analyses: function, conservation, co-regulation, dark kinases

### 2.1 How much of the phosphoproteome is functional?
- **Early scepticism.**
  - Lienhard 2008 (PMID 18603430 [V]) asked whether many sites are "non-functional".
  - Landry 2009 (PMID 19349092 [V]) showed weak evolutionary constraint on most phosphosites: many evolve nearly as fast as unmodified residues.
  - Beltrao 2012 (PMID 22817900 [V]) used about 200k PTMs from 11 species. It found PTM "hotspots" in domain families (validated on HSP70) and concluded that only a fraction of sites are likely to matter biologically [V abstract].
- **Fast turnover with functional ancient sites.** Studer 2016 (PMID 27738172 [V]) profiled 18 fungal phosphoproteomes:
  - Phosphosites diverge rapidly.
  - Ancient sites are enriched at protein interfaces and are more often functional.
  - Kinase motif frequencies and activities shifted around the whole-genome duplication.
- **Conserved hotspots.** Strumillo 2019 (PMID 31036831 [V]) found conserved phosphorylation hotspots in eukaryotic domain families.
- **Kinase specificity evolution.** Bradley 2019 (PMID 31233486 [V]) found most kinase motifs date to the last eukaryotic common ancestor, and specificity is largely unchanged after recent duplications. Bradley 2024 (PMID 38861992 [V title]) found the substrate "quality" of CK2 target sites determines their function and evolution.
- **Tyrosine loss.** Tan 2009 (PMID 19589966 [V title]) reported positive selection for tyrosine loss in metazoan evolution, consistent with avoidance of off-target pY.
- **The Ochoa 2020 functional score** (PMID 31819260 [V]):
  - Built from 6,801 re-analysed experiments, giving 119,809 human sites.
  - 59 features (proteomic, structural, regulatory, evolutionary) were integrated into a single functional score.
  - Validated sites included SMARCC2 sites in neuronal differentiation [V].
  - **Caveat for us:** the score was trained on PSP-annotated functional sites and uses "regulation across conditions" as a feature. Any analysis claiming that "functional sites are more regulated" is partly circular.
- **Structural context.** Bludau 2022 (PMID 35576205 [V title]) mapped PTMs onto AlphaFold structures. Phosphosites are enriched in disordered and accessible regions; the specific structural associations are [M].
- **Stoichiometry.**
  - Most sites have low occupancy, except in mitosis, where many reach near-full occupancy (Olsen 2010, PMID 20068231 [V title]).
  - Tyr and S/T phosphorylation have distinct regulatory natures; pTyr is low-abundance and tightly controlled (Sharma 2014, PMID 25159151 [V title]).

### 2.2 Co-regulation and kinase–kinase networks
- **Kinase co-regulation:**
  - Ochoa 2016 [V]: co-regulation predicts complexes and KS pairs, and defines logical relations between kinase pairs.
  - Invergo 2020 (PMID 32437683 [V]): supervised prediction of *signed* kinase→kinase regulation from high-throughput data, specificity profiles and structure. The network is sparse but denser than textbook pathways.
  - Invergo & Beltrao 2018 (PMID 30072490 [V]) reviews network reconstruction from phosphoproteomics.
- **Site co-regulation.**
  - SELPHI (Petsalaki 2015, PMID 25948583 [V]) and SELPHI2 (Maier 2025, MCP [V in substrate review]) use correlation-based kinase-associated networks.
  - CoPheeMap/CoPheeKSA (Jiang 2025, Nat Commun 16:2766 [V in activity review]) built a co-regulation network of 26,280 sites from 1,195 CPTAC tumours and predicted 24,015 KS links.
  - **Implication:** "phosphosite co-regulation in CPTAC predicts KS" is *not* novel. Co-regulation *modules* with condition-specific activity patterns across perturbation and tumour data jointly appear unexplored.
- **Kinase interactome.** Buljan 2020 (PMID 32707033 [V]) AP-MS'd more than 300 kinases and reported more than 5,000 new interactions, with functions for PIM3 and POMK.
- **Lun 2017** (PMID 28092656 [V]) used mass cytometry: node *abundance* shapes signaling network state and dynamics. This is single-cell support for abundance–activity coupling.

### 2.3 Dark kinases and the annotation gap
- Needham 2019 [V]:
  - <5% of sites have a known kinase.
  - About 90% of annotated sites belong to the top 20% of kinases.
  - More than 150 kinases have no known substrate.
  - **The number of substrates does not correlate with disease relevance** (human mutation and mouse phenotype).
- The Dark Kinase Knowledgebase (Berginski 2021, PMID 33079988 [V]) and the Clinical Kinase Index for prioritising understudied kinases in cancer (Essegian 2020, PMID 33205077 [V]) define the IDG "dark" set of about 160 kinases [M count].
- **Opportunity.** The Kinase Library, the Sugiyama 2019 in-vitro screen (PMID 31324866 [V]; 385 kinases, 175k pairs [V in activity review]) and our compendium together allow *activity profiles for dark kinases across ~1,900 conditions*. This has not been published as an atlas. Ochoa 2016 covered 215 mostly well-studied kinases.

### 2.4 Mutations at phosphosites and network-attacking mutations
- Reimand 2013 (PMID 24089029 [V]) mapped the mutational landscape of phosphorylation signaling in cancer (the ActiveDriver line). ActiveDriverDB covers disease and germline variants at PTM sites (Krassowski 2018, PMID 29126202 [V]).
- MIMP predicts mutations that rewire KS motifs (Wagih 2015, PMID 25938373 [V]).
- Creixell 2015 (PMID 26388441 [V]) decoded "network-attacking mutations" kinome-wide.
- Lundby 2019 (PMID 31585087 [V]) showed that mutations near pY sites switch SH2 recruitment and rewire EGFR signaling in lung tissue.

---

## 3. CPTAC phosphoproteomics: what kinase-level claims exist

CPTAC data are TMT-based, so ratio compression applies, with about 30–40k sites per cohort [M]. The pan-cancer harmonised resource covers >1,000 tumours in 10 cohorts (Li Y 2023 Cancer Cell, PMID 37582339 [V]), and LinkedOmicsKB precomputes associations for 125,969 phosphosites (Liao 2023, PMID 37619559 [V]).

**Common analysis toolkit:** PTM-SEA (Broad), KSEA/KEA3/Kinase Library (pan-cancer), BlackSheep outlier analysis for "hyperphosphorylated kinases", CausalPath (Babur 2021, PMID 34179843 [V]), and cis/trans effects of drivers.

### 3.1 Per-cohort flagship claims (verified from abstracts or full text where marked)

| Cohort | Paper | Kinase/phospho claims |
|---|---|---|
| BRCA (2016) | Mertins 2016 Nature, PMID 27251275 [V] | Besides ERBB2, amplicon-associated hyperphosphorylated kinases: CDK12, PAK1, PTK2, RIPK2, TLK2. 5q loss (CETN3, SKP1) linked to EGFR; SKP1 loss linked to SRC. A GPCR phospho-cluster not seen at mRNA level. |
| BRCA (2020) | Krug 2020 Cell, PMID 33212010 [V] | PTM-SEA CDK4/CDK6 scores correlate with proliferation (MGPS) in HR+/ERBB2− but not TNBC. Rb phosphorylation status for CDK4/6i response. Inferred mTOR activity higher in Rb-high TNBC. "Associations between tumour-suppressor loss and targetable kinases" (BlackSheep outliers, A-loop enrichment). |
| OV (2016) | Zhang H 2016 Cell, PMID 27372738 [V] | CNA effects on proteome; pathways where rearrangements converge; acetylation and HRD. Kinase specifics [M]. |
| OV (2020) | McDermott 2020 Cell Rep Med, PMID 32529193 [V] | Mitotic kinases and replication stress implicated in chromosomal instability. |
| COAD | Vasaikar 2019 Cell, PMID 31031003 [V] | Rb phosphorylation linked to proliferation and reduced apoptosis, a rationale for targeting Rb phosphorylation (CDK). Glycolysis linked to low CD8 in MSI-H. |
| ccRCC | Clark 2019 Cell, PMID 31675502 [V title; claims via full text] | Tumour vs NAT: CDK1 and MAPK1 substrates top-ranked. CDK7→MCM2 and WEE1→CDK1 elevated. CDK1 Y15 higher in aggressive tumours. EGFR up in nearly all tumours. EIF4EBP1 is a MAPK/mTOR convergence point, suggesting combined mTOR+MEK. STAT1 phosphorylation in the CD8-inflamed subtype. |
| UCEC (2020) | Dou 2020 Cell, PMID 32059776 [V full text] | CDK4 activity higher in CNV-high endometrioid tumours. CDK12 hyperphosphorylated (PARPi sensitisation). DNMT1 hyperphosphorylated in DDR-high tumours. |
| UCEC (2023) | Dou 2023 Cancer Cell, PMID 37567170 [V] | PIK3R1 in-frame indels → higher AKT phosphorylation and AKTi sensitivity. CTNNB1 hotspots cluster near S45 degron phosphosites. |
| LUAD | Gillette 2020 Cell, PMID 32649874 [V full text] | ALK fusions → outlier ALK pY1507. EGFR-mutant → PTPN11 pY62 (RNA and protein unchanged). ALK-fusion → PTPN11 pY546/pY584. KRAS-mut → SOS1 pS1161, MAPK1 and TAK1 (MAP3K7) activation. EGFR-mut → ROCK1. Kinase outliers: PRKCD (KRAS), BRAF (TP53), WEE1 (EML4-ALK). STK11 → SMAD4 pS138 and immune-cold. |
| LSCC | Satpathy 2021 Cell, PMID 34358469 [V full text] | NSD3 as alternative driver in FGFR1-amplified tumours. CDK5 protein and PTM-SEA activity up in the Classical (NRF2) subtype. RB1 phosphorylation should inform CDK4/6i. |
| HNSCC | Huang C 2021 Cancer Cell, PMID 33417831 [V full text] | Two modes of EGFR activation, with ligand abundance stratifying anti-EGFR mAb response. CHEK2 the only kinase with a hyperphosphorylated activating site. CIN subtype has high CDK4/6 activity (Rb hyperphosphorylation). |
| GBM | Wang LB 2021 Cancer Cell, PMID 33577785 [V] | PTPN11 and PLCG1 phosphorylation as RTK-pathway switches. Druggable KS outlier pairs for GSK3B, AKT1, MAPK1/3, EGFR. Follow-up Liu 2024 (PMID 38981438 [V]): PTPN11 is central across high-grade gliomas. |
| PDAC | Cao 2021 Cell, PMID 34534465 [V abstract] | Resource paper; specific kinase claims not extracted [M]. KRAS→MAPK signaling in tumour vs NAT is expected and not novel. |

### 3.2 Pan-cancer claims
- **Geffen 2023 Cell** (PMID 37582358 [V full text]): 1,110 patients, 11 cancer types. Methods were PTM-SEA, Kinase Library, CausalPath and CLUMPS-PTM. Its 33 multi-omic clusters [V via search summary] include:
  - **An HRD/acute-hypoxia cluster** with raised ATM, CHEK1 and CHEK2 PTM-SEA activity and CDK1/2/4/6 activity. Site evidence: PARP1 S782 (PKA) and S179 (ATR), EXO1 S714 (ATM), and CDK2→WRN S1133.
  - **MMR-deficient tumours** with ATM microsatellite indels and *reduced* ATM-mediated phosphorylation of PRKDC S3205 and PNKP S114/T118.
  - **Acetylation–phosphorylation crosstalk:** acetylated lysines at +3 abolish CDK1 phosphorylation (RSF1 S1375, which primes GSK3 and then PLK1 docking). In general, proximal acetylation modulates S/T kinase specificity.
  - Histone PTM regulation, including CSNK2A1-mediated sites.
- **Li Y 2023 Cell "Pan-cancer proteogenomics connects oncogenic drivers to functional states"** (PMID 37582357 [V abstract]):
  - Cis- and trans-effects of drivers at RNA, protein and phospho levels.
  - Mutations and CNAs rewire protein interaction networks.
  - **"Most cancer genes converge toward similar molecular states denoted by sequence-based kinase activity profiles."**
  - So "different drivers converge on shared kinase states" is a published claim.
- **Petralia 2024 Cell** (PMID 38359819 [V full text]): 7 pan-cancer immune subtypes. Kinase activity came from KEA3 and the Kinase Library, plus deconvolution into tumour vs immune/stromal compartments. Findings:
  - FAK (PTK2), CDK19 and CDK20 are more active in *hot* tumour cells; FAK is significant in LUAD and LSCC.
  - FYN, LYN and LCK activity is lower in cold tumours.
  - CDK1/2 are active in CD8−/IFNG± subtypes.
  - ATR activity tracks its own CNV in HNSCC.
  - **Immune-subtype × kinase activity is therefore taken.**
- **Savage 2024 Cell** (PMID 38917788 [V]): 2,863 druggable proteins. CPTAC protein overexpression or hyperactivation was combined with DepMap to identify dependencies. Also proteogenomic synthetic lethality for tumour-suppressor loss, and KRAS neoantigens. Portal: targets.linkedomics.org. **"CPTAC kinase hyperphosphorylation + DepMap dependency" is taken at the general level.**
- **CoPheeKSA** (Jiang 2025) uses CPTAC co-phosphorylation for KS prediction (see §2.2).
- **ProKAP** (Calinawan et al., bioRxiv 2021 [V via search]) is a KEA-based pan-cancer kinase activity portal. It found several kinases highly active in CD8-enriched ccRCC.
- **benchmarKIN tumour benchmark** [V]: kinase protein top/bottom 5% and activating sites were used as truth in 7 cohorts. COAD, OV and PDAC were excluded for low site–protein correlation.

### 3.3 Synthesis: what is *not* novel in CPTAC kinase analyses
- Driver-kinase cis effects (EGFR, ERBB2, ALK, PIK3R1→AKT).
- Proliferation-linked CDK activity and Rb phosphorylation.
- DNA-damage kinase activity in HRD tumours.
- PTPN11/SHP2 phosphorylation downstream of RTKs (LUAD, GBM).
- Immune-subtype kinase activities (FAK, LCK, LYN).
- Kinase outliers matched to DepMap dependencies.
- Convergence of drivers on shared kinase states.
- Co-phosphorylation-based KS prediction.

### 3.4 What remains open in CPTAC
1. Kinase *network* states (pairs and triads of kinases) rather than single-kinase activities. Examples: which kinase–kinase dependencies seen in perturbation data hold in tumours, and which break.
2. Separating abundance-driven from regulation-driven site variation, kinome-wide.
3. Using *perturbation-derived* footprints, rather than curated KS sets, to score tumours: a phospho "connectivity map".
4. Dark-kinase activity in tumours.
5. Testing whether site-level modules are conserved across cohorts and tissue-of-origin (using GTEx baselines).
6. Systematic paradox/feedback signatures in tumours treated in vivo. CPTAC is mostly treatment-naive, so this is limited.

---

## 4. Feedback, rewiring, resistance and off-target discovery

### 4.1 Classic feedback axes (mostly established by targeted assays and confirmed by phosphoproteomics)
- **MEK inhibition in TNBC:** ERK loss → c-MYC degradation → induction of multiple RTKs. Prevented by blocking MYC degradation (Duncan 2012 Cell, PMID 22500798 [V]; MIB/MS kinome profiling).
- **Lapatinib in HER2+ breast cancer:** kinome reprogramming blocked by BET inhibitors (Stuhlmiller 2015, PMID 25865888 [V]). Phosphorylation-mediated glycolysis reprogramming in lapatinib resistance (Ruprecht 2017, PMID 28209619 [V]).
- **AKT inhibition** relieves FOXO-dependent suppression of RTK expression and activity (Chandarlapaty 2011, PMID 21215704 [V]; FOXO mechanism [M]). mTOR kinase inhibition causes biphasic, feedback-dependent AKT regulation (Rodrik-Outmezguine 2011, PMID 22140653 [V]). PI3Kβ inhibition relieves feedback on PI3Kα in PTEN-mutant tumours (Schwartz 2015, PMID 25544636 [V]).
- **RAF inhibitors** relieve ERK-dependent feedback in BRAF-V600E cells, which attenuates their effect (Lito 2012, PMID 23153539 [V]). They also paradoxically transactivate RAF dimers in BRAF-WT cells (Poulikakos 2010; Hatzivassiliou 2010 [V]).
- **Growth factors** confer kinase-inhibitor resistance through phosphoproteome-level rewiring (Koch 2016, PMID 27794612 [V title]).

### 4.2 Proteome-scale inhibitor studies
- **Wilkes 2015 PNAS** (PMID 26060313 [V]):
  - Empirical PI3K/AKT/mTOR and MEK/ERK network topology from phosphoproteomics.
  - Cells chronically resistant to PI3K or mTORC1/2 inhibitors remodel differently.
  - **Cells resistant to the same drug had heterogeneous network states**, and these predicted combination responses.
- **Hijazi 2020** (PMID 31959955 [V]): kinase network topology from 61 inhibitors (§1.2).
- **Coppé 2019 NCB** (PMID 31160710 [V]): a 228-peptide kinase-activity assay covering >60 kinases. It found intrinsic BRAF-V600E resistance in CRC through parallel PDPK1 and PRKCA activation, and RPS6KB1/PIM1 as poor-outcome vulnerabilities in melanoma.
- **decryptM** (Zecha 2023 Science, PMID 36926954 [V]):
  - 31 drugs × 13 lines, 1.8M dose–response curves.
  - Drug-specific PTM signatures of kinase inhibitors.
  - DNA-damage signatures for chemotherapies.
  - Rituximab kills CD20+ B cells by *over-activating* BCR signaling.
  - Curves are classified by CurveCurator (Bayer 2023, PMID 38036588 [V]).
  - **Our compendium contains 790 decryptM contrasts.** The original paper reported curves and examples. It did not systematically compare site EC50 with target Kd_app, so it did not derive network depth or amplification at proteome scale [M; check the supplement before claiming].
- **Kinobeads target landscapes:**
  - Klaeger 2017 (PMID 29191878 [V]): 243 clinical drugs.
  - Reinecke 2024 (PMID 37904048 [V]): 1,183 compounds and 500k compound–target interactions. It found reasonably selective compounds for 72 kinases, and used phosphoproteomics to identify PKN3 substrates for engagement.
- **Tanzer 2021** (PMID 34663829 [V title]): phosphoproteomics uncovered a key role for CDKs in TNF signaling. This study is in our compendium (U937 TNF time course).
- **VESPA** (Rosenberger 2023 bioRxiv, PMID 36824919 [V]): 7-point phospho time series of CRC cells on clinical inhibitors. Kinase/phosphatase activity inference revealed adaptive rewiring, validated by CRISPR-KO.
- **Paradoxical activation-loop hyperphosphorylation by ATP-competitive inhibitors:** AKT (Okuzumi 2009 [V]); RAF (above). There is no proteome-wide catalogue of which kinases show inhibitor-induced A-loop hyperphosphorylation [M: none found].

### 4.3 Open problems
1. A **systematic, cross-context feedback atlas**: for each targeted node, which upstream sites rise (relief of negative feedback), how consistently across cell lines, and at what dose relative to Kd_app. Hill 2017 (PMID 28017544 [V]) showed causal edges are context-specific using RPPA (about 35 phosphoproteins, 32 contexts [M]). No equivalent exists at the proteome scale.
2. **Off-target discovery from phospho-footprints.** This means sites that respond at doses matching a non-annotated kinase's Kd_app (Klaeger) rather than the intended target, which flags functional off-targets.
3. **Separating direct from indirect effects** using dose and time (EC50 ordering, time-to-response).
4. **Distinguishing resistance mechanisms** that pre-exist (baseline tumour heterogeneity, CPTAC) from those induced (perturbation).

---

## 5. Temporal dynamics and causal/network inference

### 5.1 What time courses taught us
- **EGF dynamics:**
  - Olsen 2006 Cell (PMID 17081983 [V title]): the first global, site-specific EGF time course in HeLa. It showed phosphorylation dynamics are largely site-specific, even on the same protein [M].
  - Reddy 2016 (PMID 26929352 [V title]): early (seconds-scale) EGFR signaling dynamics.
  - Francavilla 2016 (PMID 27136326 [V]): ligand-dependent EGFR trafficking switches.
- **Insulin/AKT:**
  - Humphrey 2013 (PMID 23684622 [V]): AKT directly regulates mTORC2 (Sin1 phosphorylation [M]).
  - Humphrey 2015 (PMID 26280412 [V]): high-throughput in vivo insulin signaling dynamics. Mouse, not in our compendium.
- **Specific vs promiscuous phosphorylation** is resolved temporally: early responses are specific and later ones broader (Kanshin 2015, PMID 25704821 [V title]; yeast).
- **Mitosis:** Dephoure 2008 (PMID 18669648 [V]); Olsen 2010 [V]; Aurora/PLK substrate modules (Kettenbach 2011, PMID 21712546 [V]).
- **Cytokines:** CDKs are part of TNF signaling (Tanzer 2021 [V]).

### 5.2 Network and causal inference approaches
- **Prior-knowledge-constrained optimisation:**
  - PHONEMeS (Terfve 2015, PMID 26354681 [V]): signal propagation from discovery phosphoproteomics with inhibitor data.
  - PHOTON (Rudolph 2016, PMID 28009266 [V]).
  - CARNIVAL (Liu 2019, PMID 31728204 [V]).
  - COSMOS (Dugourd 2021, PMID 33502086 [V]).
  - CausalPath (Babur 2021, PMID 34179843 [V]).
- **Time-series key-event detection:** ClueR (Yang 2015, PMID 26252020 [V]); PhosR [V].
- **Data-only causal discovery:**
  - HPN-DREAM (Hill 2016, PMID 26901648 [V]): a community assessment on RPPA time courses with inhibitors. Performance was modest, and the best approaches mixed prior knowledge with data.
  - Context specificity of causal edges (Hill 2017 [V]).
  - Single-cell Bayesian networks (Sachs 2005, PMID 15845847 [V]) remain the canonical demonstration that interventional data are needed.
- **Modular response analysis** (Kholodenko 2002, PMID 12242336 [V]) infers local response coefficients from steady-state responses to perturbations of each module. It is rarely applied at phosphoproteome scale; we found no example [M].
- **EGF benchmark** (Garrido-Rodriguez 2026 [V]): literature networks outperform data-derived KS networks for pathway recovery. This is a sobering baseline for any "we reconstructed the network" claim.

### 5.3 Gaps
1. Few human time courses reach seconds-to-hours with more than 8 timepoints. Our compendium has scattered time series (EGF, TNF, FGF). Cross-study *temporal ordering* (which sites respond first, robustly across studies) has not been compiled.
2. Ultrasensitivity and amplification are classic theory (Huang & Ferrell 1996, PMID 8816754 [V]) but were never tested proteome-wide with dose series. decryptM allows it.
3. Methods to infer causal direction from *dose* (instead of time) are underdeveloped.

---

## 6. Ideas to borrow from other fields (concrete applications)

| Field / method | Original use | Concrete application here | Why it might work / risk |
|---|---|---|---|
| **ICA / iModulons** (Sastry 2019, PMID 31797920 [V]); module-detection benchmark (Saelens 2018, PMID 29545622 [V]: decomposition methods beat clustering) | Independently regulated gene modules from expression compendia | ICA on the compendium site × condition matrix (with missingness handling), giving "phosphomodules". Then annotate modules by kinase motif (Kinase Library), PSP, and perturbation drivers (which drugs/kinases activate them). Project modules onto CPTAC tumours. | Captures kinase-independent and multi-kinase programmes; does not need KS annotation. Risk: 4.6% density, so restrict to well-covered sites/conditions or use PCA with missing values. |
| **MOFA / multi-omic factor analysis** (Argelaguet 2018, PMID 29925568 [V]) | Shared factors across omics | CPTAC: factors shared by RNA, protein and phospho vs phospho-specific factors. Phospho-only factors are candidates for post-translational signaling states invisible to transcriptomics. | Shows what phospho adds beyond RNA/protein, a question raised but not resolved in the CPTAC papers. |
| **Mendelian randomisation / instrumental variables** (Davey Smith & Ebrahim 2003, PMID 12689998 [V]) | Causal effect of exposure using genetic instruments | (a) *Genetic IV*: in CPTAC, kinase cis-CNV/truncating mutation → kinase protein (first stage) → substrate phosphorylation (second stage). This gives causal kinase→site effect sizes, separates direct substrates from correlates, and tests PSP edges. (b) *Pharmacological IV*: drug dose × Kd_app (Klaeger) is an instrument for target occupancy; site response regressed on predicted occupancy across drugs sharing the target. | Natural-experiment logic; Petralia used cis-CNV→activity only descriptively. Yeast phosphoQTLs: variants affecting phosphorylation differ from those affecting abundance (Grossbach 2022, PMID 35574625 [V]). Risk: pleiotropy of CNV segments (neighbouring genes); use the arm-level CNV as a negative control. |
| **Modular response analysis / local response coefficients** (Kholodenko 2002 [V]) | Network wiring from perturbation steady states | decryptM/Hijazi: treat kinase-activity modules as nodes. Infer a signed local Jacobian from multiple inhibitor perturbations per cell line, then compare across cell lines (context-specific wiring). | Uses exactly the multi-perturbation design available. Risk: needs ≥n perturbations per n modules per context; Hijazi has 61 per line. |
| **Dose–response pharmacology / ultrasensitivity** (Huang & Ferrell 1996 [V]) | Signal amplification in cascades | For each site in decryptM, compare EC50 to the target's Kd_app. EC50 ≪ Kd means amplification or ultrasensitive coupling; EC50 ≈ Kd means a likely direct substrate or proximal site; EC50 ≫ Kd means buffering or off-target. Hill slopes give cooperativity. Test whether "depth" from target correlates with EC50 shift and Hill slope. | Quantitative, testable, and not done proteome-wide [M]. |
| **Connectivity Map** (Subramanian 2017, PMID 29195078 [V]; P100 targeted phospho: Abelin 2016, PMID 26912667 [V]) | Match disease signatures to perturbation signatures | Build data-driven perturbation footprints from the compendium (e.g., "MEKi signature", "mTORi signature") and score CPTAC tumours. This finds tumours whose phospho state mimics or opposes inhibition, giving candidate intrinsic resistance or sensitivity. | P100 covers ~96 sites only, so a proteome-scale phospho-CMap applied to tumours appears new. |
| **Footprint methods à la PROGENy** (Schubert 2018, PMID 29295995 [V]) | Pathway activity from consensus perturbation-response genes | "Phospho-PROGENy": consensus responsive sites per pathway across many studies, robust to cell type. Compare with PSP KS sets for kinase-activity accuracy (benchmarKIN). | PROGENy beat pathway-membership sets for transcriptomics; same logic may beat KS sets. |
| **Co-essentiality** (Pan 2018, PMID 29778836 [V]; Wainberg 2021, PMID 33859415 [V]; DepMap: Tsherniak 2017, PMID 28753430 [V]) | Functional modules from correlated fitness effects | Compare kinase–kinase co-regulation (compendium/CPTAC) with kinase–kinase co-essentiality (DepMap). Agreement means functional coupling; disagreement (co-regulated but not co-essential) means redundancy. This could find dark kinases co-essential with well-studied ones whose substrates co-vary. | Cheap; each layer alone is well used, the combination for kinases less so. |
| **Network control theory** (Liu 2011 Nature, PMID 21562557 [V]; human PPI controllability: Vinayagam 2016, PMID 27091990 [V]) | Driver nodes needed to steer a network | Compute driver/indispensable kinases in the signed kinase–kinase network (Invergo/Hijazi/compendium-inferred). Test whether indispensable kinases are (i) more frequently regulated (Ochoa 2016 found central kinases are), (ii) more essential (DepMap), and (iii) more often drug targets. | The structural-control claim is testable with independent data; risk of circularity with degree. |
| **Systemic-risk / shock propagation** (DebtRank: Bardoscia 2015, PMID 26091013 [V]); economic input–output analysis | Contagion of shocks through weighted networks | Score each kinase's "systemic impact": the fraction of the phosphoproteome changed when it is inhibited, from Hijazi/decryptM footprints, normalised by occupancy. Compare to network-propagation predictions; outliers are kinases with unexpected reach. | A cleanly defined empirical quantity, "kinase reach", not published as such [M]. |
| **Ecology: nestedness / specialist–generalist structure** (Bascompte 2003, PMID 12881488 [V]) | Mutualistic network architecture | Is the KS bipartite network (or the empirical drug→site network) nested, i.e., do specialist kinases' substrates form subsets of generalists'? Also: are sites targeted by many kinases more robust (fewer change upon single inhibition)? This is the redundancy/robustness hypothesis. | Robustness-via-redundancy can be tested directly with inhibitor data. |
| **Single-cell trajectories / pseudotime** (Trapnell 2014, PMID 24658644 [V]) | Ordering cells along a process | Order conditions along dose (decryptM) or time to define "pseudo-dose" trajectories. Estimate site ordering (which sites switch first) and identify branch points where responses diverge between cell lines. A "phospho-velocity" analogue: fast phospho changes vs slower protein changes in decryptM (which measured both [M]). | Dose plays the role of time with known ordering; supervised trajectories. |
| **Pareto / archetype analysis** (Hart 2015, PMID 25622107 [V]; Shoval 2012, PMID 22539553 [V]) | Trade-offs among tasks | CPTAC kinase-activity space: do tumours lie in a polytope whose vertices (archetypes) correspond to signaling "tasks" (proliferation CDK, stress/DDR, RTK–MAPK, immune)? | Gives an interpretable alternative to clustering; tests whether trade-offs constrain signaling states. |
| **Evolutionary analyses** (Studer 2016; Bradley 2019 [V]) | Site age, motif evolution | Relate site age/conservation to (i) responsiveness breadth in the compendium, (ii) dose sensitivity, (iii) CPTAC variance. Hypothesis: ancient sites are responsive in fewer, more specific conditions; young sites are noisy and abundance-driven. | Independent of PSP annotation; avoid using funscore, which includes regulation. |
| **Nested effects models** (Markowetz 2005, PMID 16159925 [V]) | Infer pathway hierarchy from subset relations of downstream effects | Inhibitor A's affected-site set ⊂ inhibitor B's → B upstream of or broader than A. This gives a hierarchy of kinases from footprint containment. | Uses footprint nesting; robust to the unknown KS map. |
| **Convergent cross mapping / Granger** (Sugihara 2012, PMID 22997134 [V]) | Causality in short ecological time series | Apply to the few dense phospho time courses to test the directionality of kinase-activity coupling. | Likely underpowered (≤10 timepoints); low priority. |

---

## 7. Open questions and gaps addressable in days with our data

Data: compendium of 1,937 contrasts (decryptM 790, Hijazi 183, many time courses); CPTAC 10 cohorts; DepMap; GTEx; interactomes; Kinase Library; Ochoa functional scores; PSP; Klaeger.

1. **Proteome-wide dose-sensitivity map.**
   - Method: for each decryptM site, compute EC50 relative to the drug's on-target Kd_app.
   - Question: does EC50 shift with network distance? Are there "amplifier" sites with EC50 ≪ Kd?
   - Novelty: decryptM reported curves, but ultrasensitivity and amplification have not been quantified proteome-wide against Kd [M].
2. **Functional off-target detection.**
   - Method: sites whose EC50 matches the Kd_app of a non-intended kinase rather than the intended target.
   - Novelty: off-targets are known from binding assays (Klaeger/Reinecke); here they would be validated *functionally in cells*.
3. **Cross-context feedback atlas.**
   - Method: for each inhibited kinase family, sites whose phosphorylation *rises*, how reproducibly across 13 decryptM lines and 3 Hijazi lines, and at what dose relative to Kd.
   - Novelty: known feedbacks are axis-specific; no catalogue of the cell-line invariance of feedback exists.
4. **Paradoxical target hyperphosphorylation catalogue.**
   - Question: which kinases' activation-loop or hydrophobic-motif sites increase when the kinase is inhibited?
   - Novelty: known only for AKT and RAF case studies.
5. **Phosphomodules by ICA on the compendium, projected onto CPTAC.**
   - Novelty: there are kinase-level atlases (Ochoa 2016) and CPTAC co-phosphorylation networks (CoPheeMap), but no perturbation-anchored site modules transferred to tumours.
6. **Phospho-CMap.**
   - Method: score CPTAC tumours with data-driven inhibitor footprints and find tumours "mimicking" inhibition (e.g., low-MEK-footprint, KRAS-mutant tumours).
   - Novelty: P100-CMap used ~96 sites in cell lines only.
7. **Genetic-IV (MR-like) kinase→substrate effects in CPTAC.**
   - Method: cis-CNV and truncating mutations of kinases as instruments.
   - Novelty: gives quantitative in-vivo causal KS evidence. Petralia used cis-CNV descriptively, for ATR only.
8. **Which perturbation-derived kinase–kinase relations hold in tumours?**
   - Method: test Hijazi/compendium kinase–kinase edges as co-variation in CPTAC, after abundance adjustment.
   - Novelty: rewiring between cell lines and tumours has not been quantified; Hijazi applied it only to AML and breast.
9. **Abundance- vs regulation-driven site variance, kinome-wide.**
   - Method: in CPTAC, partition each site's variance into protein-abundance, kinase-activity and residual components. Ask which kinases' substrates are dominated by stoichiometry changes.
   - Novelty: benchmarKIN hinted that protein normalisation hurts, but there is no systematic decomposition.
10. **Dark-kinase activity atlas.**
    - Method: Kinase Library motifs plus perturbation co-regulation, giving activity profiles for dark kinases across ~1,900 conditions and CPTAC. Validate against DepMap co-essentiality and GTEx tissue expression.
    - Novelty: Ochoa 2016 covered 215 mostly well-studied kinases.
11. **Site age vs regulatory behaviour.**
    - Question: do evolutionarily ancient sites show more context-specific and dose-sensitive responses than young sites?
    - Novelty: Studer showed ancient sites are more functional; their *dynamic behaviour* in human perturbations is untested.
12. **Kinase "reach" (systemic impact).**
    - Method: the fraction of the phosphoproteome perturbed per unit of target occupancy. Compare with network centrality, controllability and essentiality.
    - Novelty: a new empirical quantity linking network theory to data.
13. **Footprint nesting hierarchy (NEM-style) from inhibitor data.**
    - Method: infer kinase ordering from containment of affected-site sets.
    - Novelty: NEMs have been used for RNAi transcriptomics, not phospho inhibitor panels [M].
14. **Kinase co-regulation vs co-essentiality.**
    - Method: compare kinase–kinase co-regulation (compendium/CPTAC) with DepMap co-essentiality. Discordance flags redundancy or compensation.
    - Novelty: each layer has been studied separately [M].
15. **Pan-cancer archetypes of kinase activity.**
    - Method: Pareto analysis of CPTAC kinase activities.
    - Novelty: clustering-based subtypes exist (Geffen, Petralia); trade-off geometry has not been tested.
16. **Phosphatase footprints.**
    - Method: infer PP2A-B56 (LxxIxE docking), PP1 and calcineurin activity from motif-defined substrate sets in perturbations and tumours.
    - Novelty: phosphatases are absent from activity atlases.
17. **Temporal ordering across studies.**
    - Method: meta-analyse the compendium's time courses to define "early robust responder" sites for EGF, TNF and growth factors. Test whether early sites are enriched for direct KS motifs (Kinase Library), following Kanshin's specific-early, promiscuous-late idea.
    - Novelty: Kanshin's idea has been shown in yeast, not tested in human meta-analysis.
18. **Tissue-of-origin vs tumour signaling.**
    - Method: GTEx kinase expression baselines vs CPTAC activities. Are tumour kinase activities explained by lineage? Find lineage-independent activity programmes.
    - Novelty: CPTAC papers compare against NAT per cohort but rarely against the normal-tissue kinome expression baseline across cohorts.
19. **Interactome-supported activity inference for dark kinases.**
    - Method: Buljan AP-MS, BioPlex or OpenCell partners × co-regulated sites.
    - Novelty: tests whether physical partners are better activity readouts than motifs.
20. **Acetyl–phospho crosstalk generalisation.**
    - Question: does +3K or other lysine context at CDK/basophilic motifs predict reduced phosphorylation responsiveness in the compendium (where acetylation is not measured)?
    - Novelty: extends Geffen's single-motif biochemistry. Weak with our data, so low priority.

### Pitfalls to remember when claiming novelty
- **Circularity:**
  - Ochoa funscore uses regulation features.
  - CoPheeKSA and SELPHI2 are trained on CPTAC or compendium co-regulation.
  - Hijazi PDTs are derived from the Hijazi data.
  - PSP in-vivo edges may come from compendium studies.
- **Paralogs:** motif-based "activity" is family-level. Claim kinase-specific results only with orthogonal support (protein abundance, A-loop site, DepMap).
- **CPTAC sample composition:** tumour purity and immune infiltration drive many kinase activities (SYK, LCK, BTK). Adjust for purity, ESTIMATE scores or deconvolution before claiming tumour-intrinsic signaling.
- **TMT ratio compression** and cohort batch effects call for per-cohort analysis with meta-analysis across cohorts.

---

## References (verified PMIDs; see tags in text)
- Basanta 2023 PMID 37352361
- Bardoscia 2015 PMID 26091013
- Bascompte 2003 PMID 12881488
- Beltrao 2012 PMID 22817900
- Babur 2021 PMID 34179843
- Badia-i-Mompel 2022 PMID 36699385
- Bayer 2023 PMID 38036588
- Berginski 2021 PMID 33079988
- Bludau 2022 PMID 35576205
- Bradley 2019 PMID 31233486
- Bradley 2024 PMID 38861992
- Buljan 2020 PMID 32707033
- Cao 2021 PMID 34534465
- Casado 2013 PMID 23532336
- Chandarlapaty 2011 PMID 21215704
- Clark 2019 PMID 31675502
- Coppé 2019 PMID 31160710
- Creixell 2015 PMID 26388441
- Crowl 2022 PMID 35879309
- Davey Smith 2003 PMID 12689998
- Dephoure 2008 PMID 18669648
- Dou 2020 PMID 32059776
- Dou 2023 PMID 37567170
- Dugourd 2021 PMID 33502086
- Duncan 2012 PMID 22500798
- Essegian 2020 PMID 33205077
- Francavilla 2016 PMID 27136326
- Frejno 2020 PMID 32686665
- Garrido-Rodriguez 2026 PMID 41741430
- Geffen 2023 PMID 37582358
- Gillette 2020 PMID 32649874
- Grossbach 2022 PMID 35574625
- Hart 2015 PMID 25622107
- Hatzivassiliou 2010 PMID 20130576
- Hernandez-Armenta 2017 PMID 28200105
- Hijazi 2020 PMID 31959955
- Hill 2016 PMID 26901648
- Hill 2017 PMID 28017544
- Huang C 2021 PMID 33417831
- Huang & Ferrell 1996 PMID 8816754
- Humphrey 2013 PMID 23684622
- Humphrey 2015 PMID 26280412
- Invergo 2018 PMID 30072490
- Invergo 2020 PMID 32437683
- Johnson 2023 PMID 36631611
- Kanshin 2015 PMID 25704821
- Kettenbach 2011 PMID 21712546
- Kholodenko 2002 PMID 12242336
- Kim 2021 PMID 33626354
- Klaeger 2017 PMID 29191878
- Koch 2016 PMID 27794612
- Koenig 2026 PMID 42816593
- Kong 2025 PMID 39838083
- Krassowski 2018 PMID 29126202
- Krug 2019 PMID 30563849
- Krug 2020 PMID 33212010
- Landry 2009 PMID 19349092
- Li Y 2023 Cancer Cell PMID 37582339
- Li Y 2023 Cell PMID 37582357
- Liao 2023 PMID 37619559
- Lienhard 2008 PMID 18603430
- Lito 2012 PMID 23153539
- Liu A 2019 PMID 31728204
- Liu J 2024 PMID 38981438
- Liu YY 2011 PMID 21562557
- Lun 2017 PMID 28092656
- Lundby 2019 PMID 31585087
- Lussana 2024 PMID 39563468
- Markowetz 2005 PMID 16159925
- McDermott 2020 PMID 32529193
- Mertins 2016 PMID 27251275
- Mischnik 2016 PMID 26628587
- Müller-Dott 2025 PMID 40404650
- Needham 2019 PMID 30670635
- Needham 2024 PMID 39577414
- Ochoa 2016 PMID 27909043
- Ochoa 2020 PMID 31819260
- Okuzumi 2009 PMID 19465931
- Olsen 2006 PMID 17081983
- Olsen 2010 PMID 20068231
- Pan 2018 PMID 29778836
- Petralia 2024 PMID 38359819
- Petsalaki 2015 PMID 25948583
- Poulikakos 2010 PMID 20179705
- Reddy 2016 PMID 26929352
- Reimand 2013 PMID 24089029
- Reinecke 2024 PMID 37904048
- Rodrik-Outmezguine 2011 PMID 22140653
- Rosenberger 2023 PMID 36824919
- Rudolph 2016 PMID 28009266
- Ruprecht 2017 PMID 28209619
- Sachs 2005 PMID 15845847
- Saelens 2018 PMID 29545622
- Sastry 2019 PMID 31797920
- Satpathy 2021 PMID 34358469
- Savage 2024 PMID 38917788
- Schubert 2018 PMID 29295995
- Schwartz 2015 PMID 25544636
- Sharma 2014 PMID 25159151
- Shoval 2012 PMID 22539553
- Strumillo 2019 PMID 31036831
- Studer 2016 PMID 27738172
- Stuhlmiller 2015 PMID 25865888
- Subramanian 2017 PMID 29195078
- Sugihara 2012 PMID 22997134
- Sugiyama 2019 PMID 31324866
- Tan 2009 PMID 19589966
- Tanzer 2021 PMID 34663829
- Terfve 2015 PMID 26354681
- Trapnell 2014 PMID 24658644
- Tsherniak 2017 PMID 28753430
- Vasaikar 2019 PMID 31031003
- Vinayagam 2016 PMID 27091990
- Wagih 2015 PMID 25938373
- Wainberg 2021 PMID 33859415
- Wang LB 2021 PMID 33577785
- Wilkes 2015 PMID 26060313
- Wilkes 2017 PMID 28674151
- Wiredja 2017 PMID 28655153
- Yaron-Barir 2024 PMID 38720073
- Yang 2015 PMID 26252020
- Yılmaz 2021 PMID 33608514
- Zecha 2023 PMID 36926954
- Zhang H 2016 PMID 27372738

Not verified this session:
- Engholm-Keller 2019 (KinSwingR, PLoS Biol)
- Beekhof 2019 (INKA, MSB 15:e8250)
- Jiang 2025 (CoPheeKSA, Nat Commun 16:2766; verified in the activity review)
- Maier 2025 (SELPHI2, MCP)
- Calinawan ProKAP (bioRxiv 2021)
- Abelin 2016 PMID 26912667 (verified)
