"""I26 step 8: curated misattribution / correction list (top evidence-ranked, verified-inhibition alternatives,
<=2 groups supporting the site, alt mean z<=-1.5, motif>=85 for alt member or a canonical downstream kinase) + literature notes
(web check 2026-10-29; PSP KS flat file has no PMIDs, so original references were located by literature search)."""
import pandas as pd
R26 = "/workspace/results/b_I26"
R = pd.read_csv(f"{R26}/tables/corrections.tsv", sep="\t")
v = R[(R.evidence == "verified") & (R.n_groups_supporting_site <= 2) & (R.alt_mean_z <= -1.5)
      & (R[["alt_motif_pct", "downstream_motif_pct"]].max(axis=1) >= 85)].drop_duplicates(["site", "psp_group"])
LIT = {
 "P29317_S897": ("supports correction", "Original AKT site (Miao 2009, Cancer Cell). Later work: EGF- and radiation-induced S897 blocked by MEK/RSK inhibitors, NOT by PI3K/AKT inhibitors (Hamaoka 2016 Cell Signal; 2017-18); RSK in vitro (Zhou 2015 Mol Cell). PSP already lists RSK1/3 -> AKT edge context-specific."),
 "Q9UJU6_S269": ("no kinase evidence for AKT found", "Li 2014 Oncotarget: S269/T291 14-3-3 sensor sites, phosphatase-sensitive; kinase not identified in abstract/snippets. AKT1 attribution unsupported by our data; MEK/ERK-RSK pathway dependence (5/5 studies) with RSK motif 94th pct. Candidate novel RSK site."),
 "P01106_S77": ("unclear", "GSK3 canonical site is T58; S77 is not a recognised GSK3 site. PSP also lists MTOR; data support mTOR-pathway dependence."),
 "Q14247_S405": ("supports correction", "Cortactin S405/S418 are classic ERK sites (Campbell 1999); PSP PRKCD/CDK1/2 edges not supported; ERK already in PSP."),
 "O00443_S259": ("original = cellular correlation", "Didichenko 2003: S259 phosphorylated at G2/M (Cdc2-linked) and by JNK under UV stress; acute CDK inhibition in asynchronous cells does not lower it, MEK inhibition does (5/5 studies, ERK motif 90). Context: mitosis-specific vs basal ERK-dependent."),
 "Q01167_S373": ("original mapped other sites", "Marais 2010 JBC mapped CDK sites S368/S423 (S373 appears as database inference, SIGNOR low score). Our data: mTOR-inhibitor dependent (FOXK1/2 dephosphorylation under mTORC1 control is documented for FOXK1; He 2018 Mol Cell)."),
 "P40763_S727": ("multi-kinase site", "STAT3 S727: >20 PSP kinases; ERK already annotated. MTOR edge not supported (4 studies); ERK supported 3/3."),
 "P36507_T394": ("supports correction", "Vendor literature attributes T394 to CDK5 (Banks 2015) but also to ERK2 negative feedback (analogous MEK1 T292). Our data: RAF-inhibitor dependent, ERK motif 97th pct -> ERK feedback site."),
 "P27708_S1859": ("pathway-consistent", "Canonical S6K1 site (Ben-Sahra 2013, Robitaille 2013); S6K group only 1 verified study -> treat as untested, mTOR support consistent. NOT a misattribution."),
 "Q9Y613_S498": ("pathway-consistent", "EGFR-inhibitor dependent; ERK annotation consistent at pathway level (ERKi studies weak)."),
 "P35568_S1101": ("multi-kinase, context", "S6K1 in vitro + S6K1-KO liver (Tremblay 2007 PNAS); PKCtheta (Li 2004 JBC). Acute data: ERK-inhibitor dependent (3/3), RSK motif 100 -> RSK also plausible; S6K edge untestable-ish (1 verified S6Ki study)."),
 "O14757_S280": ("supports correction", "AKT attribution (Puc 2005; Shtivelman 2002) contested: Li 2012 MBoC showed p90RSK, not AKT, phosphorylates S280 stoichiometrically. Our data agree (ERKi 2/2, RSK3 motif 95). RSK1 already in PSP."),
 "Q9H467_S110": ("context (mitotic)", "Gao 2011 NCB: Cdk1 phosphorylates CUEDC2 in mitosis (S110 per vendor). Interphase data: RAF/ERK-pathway dependent (ERK motif 93). Likely both: mitotic CDK1, interphase ERK."),
 "Q9UN36_S332": ("supports correction", "NDRG2 S332: PSP lists AKT1/RSK1/S6K/SGK1/PKCtheta (in vitro panel). Data: ERK-pathway (3/3), RSK1 motif 97 -> RSK."),
 "Q9ULW0_T369": ("candidate", "MAPK1 in vitro only; data CDK-dependent (CDK5/CDK1 motif 93)."),
 "P25098_S670": ("both kinases documented", "ERK1 site in vitro/feedback (Pitcher 1999; Elorza 2000); CDK2 at G2/M (Penela 2010 PNAS). Data: ERK-dependent, not CDK2 -> CDK2 edge is mitosis/cell-cycle-context-specific."),
 "P0C1Z6_S180": ("no literature found", "No TFPT S180 kinase report found; CDK1 edge likely from in vitro screen. Data: RAF-inhibitor dependent, ERK motif 100."),
 "P35568_S629": ("candidate", "AKT site (Luo 2007). Data: MEK-dependent (2/2), MSK1/RSK motif 98."),
 "Q9UKV8_S387": ("multi-kinase", "AGO2 S387: AKT3 (Horman 2013) and p38/MK2 (Zeng 2008). Data: ALK/ROS1-inhibitor (RSK motif 96) dependent; AKT not supported."),
 "P56524_S632": ("context", "HDAC4 S632 (14-3-3 site, CaMK/PKD family). PKCalpha edge not supported (z +4.5); mTOR-inhibitor dependent (SGK1 motif 93)."),
 "P26358_S143": ("original = in vitro + KD", "Esteve 2011 NSMB: AKT1 phosphorylates S143 in vitro; AKT1 knockdown/DN in cells. Data: AKT inhibitors no effect (2 studies), MEK-dependent (RSK3 motif 98)."),
 "P04792_S82": ("known MK2 site", "HSPB1 S82 canonical MAPKAPK2 site; AKT1 edge not supported. P38/MK2 untestable here (no verified p38 inhibition)."),
 "Q96F86_S161": ("conflicting literature", "AKT (Larance 2010 MCP, insulin) vs PIM1/3 (Bearss 2021 EMBO Rep). Data: AKT inhibitors no effect, ALK/ROS1-inhibitor dependent, RSK motif 100 -> RSK/PIM-like basophilic kinase."),
 "Q9UKA4_T1100": ("no literature found", "No AKAP11 T1100 kinase report found; data ERK-dependent (ERK1 motif 88)."),
 "P29590_S36": ("candidate", "PML S36 MAPK1/3/HIPK2; data AKT-inhibitor dependent (2/2)."),
 "Q9Y2V2_S52": ("supports correction", "CARHSP1 S52: AKT1/RSK1/SGK1 in PSP; data RSK-pathway (ALK/ROS1) dependent, AKT not."),
}
v["lit_verdict"] = v.site.map(lambda s: LIT.get(s, ("not checked", ""))[0])
v["lit_note"] = v.site.map(lambda s: LIT.get(s, ("", ""))[1])
v["correction_class"] = v.apply(lambda r: "pathway-consistent (not misattribution)" if r.same_axis else
                                ("alternative already in PSP (drop PSP edge in this context)" if (r.alt_already_psp or any(k in str(r.all_psp_kinases_site).split(";") for k in [str(r.downstream_motif_kinase), str(r.alt_motif_member)]))
                                 else "novel alternative"), axis=1)
cols = ["gene", "site", "window", "psp_kinases", "evidence_psp", "LT_LIT", "psp_n_meas", "psp_mean_z", "psp_motif", "alt_group", "alt_n_dec", "alt_n_meas",
        "alt_mean_z", "alt_motif_member", "alt_motif_pct", "downstream_motif_kinase", "downstream_motif_pct", "n_groups_supporting_site",
        "correction_class", "all_psp_kinases_site", "lit_verdict", "lit_note"]
v[cols].to_csv(f"{R26}/tables/corrections_curated.tsv", sep="\t", index=False, float_format="%.2f")
print(v.correction_class.value_counts().to_dict()); print(v.lit_verdict.value_counts().to_dict())
print(v[["gene","site","psp_kinases","alt_group","correction_class","lit_verdict"]].to_string())
