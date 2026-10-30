"""I03 step 9: figures + curated top-20 table.
fig1_validation.png  macro AUROC per score on PSP-TEST (held-out proteins) and external curated non-PSP edges
fig2_conditional_OR.png  odds ratios of new evidence among high-prior pairs (+ shuffled controls)
fig3_stratified.png  PSP rate by kinaseSub tier x replicated drop (TEST proteins)
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES = "/workspace/results/b_I03"; T = f"{RES}/tables"; import os; os.makedirs(f"{RES}/figures", exist_ok=True)
BLUE, ORANGE, GREY, INK, MUTED = "#2a78d6", "#eb6834", "#9a9a94", "#1f1f1e", "#6b6a64"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.axisbelow": True})
V = pd.read_csv(f"{T}/val_metrics.tsv", sep="\t")
pick = [("Cantley motif", GREY), ("PPI", GREY), ("perturbation: replicated drop", BLUE), ("perturbation: replicated drop (2nd study)", BLUE),
        ("tumour co-variation z", BLUE), ("NEW (LR dev)", BLUE), ("NEW (PSP-trained LR)", BLUE), ("NEW+motif (LR dev)", BLUE),
        ("NEW+motif (PSP-trained LR)", BLUE), ("kinaseSub DEV-trained, no-phospho", GREY), ("kinaseSub DEV + PPI + motif (LR)", GREY),
        ("INTEGRATED NEW+motif+kinaseSubDEV+PPI (LR dev)", ORANGE)]
ROWS = [("Cantley motif", ["Cantley motif"], GREY), ("PPI", ["PPI"], GREY),
        ("perturbation: replicated drop", ["perturbation: replicated drop", "perturbation: replicated drop (2nd study)"], BLUE),
        ("tumour co-variation", ["tumour co-variation z"], BLUE), ("NEW (pert+dose+tumour, LR)", ["NEW (LR dev)", "NEW (PSP-trained LR)"], BLUE),
        ("NEW + motif (LR)", ["NEW+motif (LR dev)", "NEW+motif (PSP-trained LR)"], BLUE),
        ("kinaseSub (DEV-trained, no phospho)", ["kinaseSub DEV-trained, no-phospho"], GREY),
        ("kinaseSub-DEV + PPI + motif (LR)", ["kinaseSub DEV + PPI + motif (LR)"], GREY),
        ("INTEGRATED: all of the above (LR)", ["INTEGRATED NEW+motif+kinaseSubDEV+PPI (LR dev)"], ORANGE)]
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2))
for i, (ax, lab, ttl) in enumerate([(axs[0], "PSP TEST proteins (held out)", "PSP edges, held-out TEST proteins (18 kinases)"),
                     (axs[1], "external curated non-PSP", "External curated non-PSP edges (29 kinases)")]):
    d = V[V.label == lab].set_index("score")
    vals = [next((d.loc[a, "macro_auroc"] for a in al if a in d.index), np.nan) for _, al, _ in ROWS]
    y = np.arange(len(ROWS))
    ax.barh(y, np.array(vals) - 0.5, left=0.5, color=[c for *_, c in ROWS], height=0.6)
    for yi, v in zip(y, vals): ax.text(v + 0.005, yi, f"{v:.2f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(y, [n for n, *_ in ROWS] if i == 0 else [""] * len(ROWS)); ax.set_xlim(0.5, 0.95)
    ax.set_title(ttl, fontsize=9.5, color=INK, loc="left"); ax.set_xlabel("macro AUROC (within-kinase ranking)"); ax.invert_yaxis()
fig.text(0.01, 0.01, "blue = new evidence (perturbation / dose / tumour); grey = sequence / PPI / kinaseSub; orange = integrated", fontsize=8, color=MUTED)
plt.tight_layout(rect=(0, 0.04, 1, 1)); plt.savefig(f"{RES}/figures/fig1_validation.png", dpi=160); plt.close()

O = pd.read_csv(f"{T}/val_conditional_OR.tsv", sep="\t")
O = O[O.term.isin(["drop2", "tum_pos", "dose_pos"])]
lab = {"drop2": "decrease ≥2 studies", "tum_pos": "tumour z>3", "dose_pos": "dose r>0.3"}
short = {"PSP, all proteins": "PSP all", "PSP, TEST proteins": "PSP TEST", "external curated non-PSP": "external", "PSP, all proteins [shuffled control]": "PSP all, shuffled", "external curated non-PSP [shuffled control]": "external, shuffled"}
fig, ax = plt.subplots(figsize=(7.5, 4.6))
rows = []
for L in ["PSP, all proteins", "PSP, TEST proteins", "external curated non-PSP", "PSP, all proteins [shuffled control]", "external curated non-PSP [shuffled control]"]:
    for t in ["drop2", "tum_pos", "dose_pos"]:
        r = O[(O.label == L) & (O.term == t)]
        if len(r): rows.append((f"{lab[t]} | {short[L]}", r.OR.iloc[0], r.lo.iloc[0], r.hi.iloc[0], "shuffled" in L))
for i, (n, o, lo, hi, sh) in enumerate(rows):
    c = GREY if sh else BLUE
    ax.plot([max(lo, 0.1), min(hi, 40)], [i, i], color=c, lw=2); ax.plot(o, i, "o", color=c, ms=7, mec="white", mew=1.5)
ax.axvline(1, color=MUTED, lw=1); ax.set_xscale("log"); ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=7.5)
ax.invert_yaxis(); ax.set_xlabel("odds ratio of being a known edge (95% CI), \nadjusted for kinase, kinaseSub tier, motif, coverage")
ax.set_title("New evidence among high-prior pairs (kinaseSub-DEV pct ≥95)", fontsize=10, loc="left", color=INK)
plt.tight_layout(); plt.savefig(f"{RES}/figures/fig2_conditional_OR.png", dpi=160); plt.close()

S = pd.read_csv(f"{T}/val_stratified_psp_rate_TEST.tsv", sep="\t")
tiers = ["95-99", "99-99.7", ">99.7"]; perts = [("no drop", GREY), ("drop in 1 study", "#86b6ef"), ("drop in >=2 studies", BLUE)]
fig, ax = plt.subplots(figsize=(6, 3.6)); w = 0.26
for j, (p, c) in enumerate(perts):
    v = [S[(S.ks_tier == t) & (S.pert == p)].psp_rate.iloc[0] * 100 for t in tiers]
    n = [S[(S.ks_tier == t) & (S.pert == p)].n.iloc[0] for t in tiers]
    x = np.arange(3) + (j - 1) * w
    ax.bar(x, v, w - 0.03, color=c, label=p)
    for xi, vi, ni in zip(x, v, n): ax.text(xi, vi + 0.3, f"n={ni}", ha="center", fontsize=6.5, color=MUTED)
ax.set_xticks(range(3), [f"kinaseSub-DEV pct {t}" for t in tiers]); ax.set_ylabel("% pairs that are PSP edges")
ax.legend(frameon=False, fontsize=8); ax.set_title("Held-out TEST proteins: replicated inhibitor-induced decrease", fontsize=10, loc="left", color=INK)
plt.tight_layout(); plt.savefig(f"{RES}/figures/fig3_stratified_TEST.png", dpi=160); plt.close()

# ---- curated top-20
CUR = [
 ("MTOR", "Q13131_S356", "PRKAA1 S356 (=AMPKa1 S347 legacy numbering)", "Rediscovery (positive control): mTORC1 phosphorylates AMPKa1 S347/a2 S345 (Ling 2020 Nat Metab) but edge absent from PSP/OmniPath/iPTMnet", "mTORC1->AMPK inhibition under nutrient stress"),
 ("MTOR", "Q9HB90_S15", "RRAGC S15", "Novel site; neighbouring RagC S21 is a reported mTORC1 feedback site (Yang 2019 EMBO J)", "Rag GTPase feedback on lysosomal mTORC1"),
 ("MTOR", "Q92609_S730", "TBC1D5 S730", "Not reported; TBC1D5 known to be dephosphorylated on mTOR inhibition (site unspecified)", "Retromer/Rab7 GAP; autophagy, ATG9 trafficking"),
 ("MTOR", "Q6IAA8_S98", "LAMTOR1 S98", "Not reported (2026 thesis title links LAMTOR1 phosphorylation to mTORC1, site unknown)", "Ragulator anchor; lysosome positioning"),
 ("MAPK3", "Q9Y6I3_S435", "EPN1 S435 (ERK1/2)", "Not reported (only mitotic CDK1 phosphorylation of epsin known)", "Clathrin endocytosis of EGFR/RTKs - ERK feedback on receptor internalisation"),
 ("MAPK3", "Q5PRF9_S271", "SAMD4B S271 (ERK1/2)", "Not reported", "Smaug2 RNA-binding translational repressor"),
 ("MAPK1", "Q6PKG0_T376", "LARP1 T376 (ERK1/2)", "Not checked in depth; LARP1 known mTORC1/S6K substrate at other sites", "TOP-mRNA translation control"),
 ("AKT2", "O15231_S64", "ZNF185 S64 (AKT1-3)", "Not reported", "LIM-domain actin-associated protein; high Ochoa functional score (0.91)"),
 ("AKT2", "Q9NZN5_S1288", "ARHGEF12/LARG S1288 (AKT1-3)", "CONFLICT/extension: reported as RSK2 site (JBC 2025); AKT dependence + MK2206 dose coherence (r=0.45) suggest AKT as second AGC kinase", "RhoA activation, invasion"),
 ("CDK1", "Q86XL3_S512", "ANKLE2 S512", "Not reported", "PP2A-ANKLE2 nuclear envelope reassembly (mitotic-exit confound)"),
 ("ATM", "Q12888_S580", "TP53BP1 S580", "53BP1 is an ATM substrate (many SQ sites) but S580 not reported", "DSB repair pathway choice"),
 ("ATR", "Q14566_S13", "MCM6 S13", "Partially reported: listed as ATR site in an MCM review table (primary ref unchecked); not in DBs", "Replication-stress response at origins/forks"),
 ("ATR", "Q9HAW4_S950", "CLSPN S950", "Not reported", "Claspin, ATR-CHK1 mediator"),
 ("ATR", "Q9UKV3_S243", "ACIN1 S243", "Not reported", "Apoptotic chromatin condensation / EJC"),
 ("CSNK2A1", "O00629_S60", "KPNA4 S60 (CK2a/a')", "CONFLICT: S60 reported as MAPK-dependent (Oncogene 2019, U0126); here CK2 inhibitors drop it in 2 studies, MEK inhibitors do not", "Importin-a3 IBB domain, nuclear import (NF-kB) in HNSCC"),
 ("CSNK2A1", "P51116_S453", "FXR2 S453 (CK2a/a')", "Not reported; CK2-FXR2 interaction in BioGRID (AP-MS)", "Fragile-X related RNA-binding protein"),
 ("GSK3B", "Q32MZ4_S116", "LRRFIP1 S116 (GSK3a/b)", "Not reported", "Wnt/b-catenin modulator, nucleic-acid sensing"),
 ("ROCK1", "O14974_S903", "PPP1R12A/MYPT1 S903 (ROCK1/2)", "Not reported (ROCK sites T696/T853 known)", "Myosin phosphatase - actomyosin contractility"),
 ("CDK9", "Q9UHB7_S1043", "AFF4 S1043 (CDK family)", "Not reported (P-TEFb phosphorylation of AFF4 known, sites T220/S212)", "Super elongation complex, transcription"),
 ("IKBKB", "P43487_S32", "RANBP1 S32", "Not reported", "Ran GTPase cycle; NF-kB-linked nucleocytoplasmic transport"),
]
E = pd.read_csv(f"{T}/edges_ranked_all.tsv.gz", sep="\t").set_index(["kinase", "site"])
D = pd.read_csv(f"{T}/top_edges_evidence_detail.tsv", sep="\t")
out = []
for k, s, name, lit, bio in CUR:
    r = E.loc[(k, s)]
    dd = D[(D.kinase == k) & (D.site == s)]
    studies = "; ".join(f"PMID{int(x.study)} {x.agents} ({x.cells}) z={x.mean_z:.1f} log2FC={x.mean_log2fc:.2f}" for x in dd.itertuples())
    out.append(dict(kinase=k, site=s, edge=name, n_studies=int(r.n_st), n_studies_decrease=int(r.n_st_drop), stouffer_drop=round(r.stouffer, 2),
                    p_coverage_matched=round(r.p_cov, 3), kinase_rank_inhibitor_drop=r.rank_inh, kinaseSub_DEV_pct=round(r.ks_dev_pct, 2),
                    kinaseSub_release_pct=round(r.ks_rel_pct, 2), motif_pct=r.motif, ppi=r.ppi, tumour_z=round(r.tum_z, 2) if pd.notna(r.tum_z) else None,
                    tumour_cohorts=r.tum_ncoh, dose_r=round(r.dose_r, 2) if pd.notna(r.dose_r) else None,
                    family_best_other=r.fam_best_other, family_other_ksdev=r.fam_other_ksdev, orphan_site=r.orphan_site,
                    other_known_kinases=r.site_known_kinases, ochoa_functional=r.functional_score, posterior_prob=round(r.post_prob, 3),
                    literature=lit, biology=bio, per_study=studies))
pd.DataFrame(out).to_csv(f"{T}/top20_curated_edges.tsv", sep="\t", index=False)
print(pd.DataFrame(out)[["edge", "n_studies", "n_studies_decrease", "stouffer_drop", "kinaseSub_DEV_pct", "motif_pct", "tumour_z", "dose_r", "posterior_prob"]].to_string())
