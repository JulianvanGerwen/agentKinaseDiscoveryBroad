"""val_I03 figures: fig_T1_mtor.png, fig_T2_kpna4.png, fig_T3_epn1.png"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
T = "/workspace/results/val_I03/tables"; F = "/workspace/results/val_I03/figures"
COL = {"direct_mTORC1": "#2a78d6", "S6K": "#eb6834", "AKT": "#1baf7a", "cand": "#222222", "ERK": "#eda100", "other": "#9a9a94"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#888",
                     "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": .6})
U = pd.read_csv(f"{T}/unit_effects.tsv", sep="\t", dtype={"study_pmid": str})
M = pd.read_csv(f"{T}/matched_ratios.tsv", sep="\t", dtype={"study_pmid": str})
D = pd.read_csv(f"{T}/dose_pEC50.tsv", sep="\t")
CP = pd.read_csv(f"{T}/cptac_meta.tsv", sep="\t")

def gcol(g):
    if g.startswith("cand") or g == "posctrl_mTOR": return COL["cand"]
    if g == "neighbour": return COL["cand"]
    return COL.get(g, COL["ERK"] if g.startswith("ERK") else COL["other"])

# ---------------- T1 ----------------
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
sites = ["RRAGC S15", "LAMTOR1 S98", "LAMTOR1 S97", "TBC1D5 S730", "PRKAA1 S356",
         "4EBP1 T46", "4EBP1 S65", "PRAS40 S183", "LARP1 S766", "ULK1 S758",
         "EIF4B S422", "RPS6 S240", "PRAS40 T246", "GSK3B S9", "TSC2 S939"]
classes = ["cat_mTORi", "rapalog", "S6Ki", "AKTi", "PI3Ki", "MEKi"]
CE = pd.read_csv(f"{T}/class_effects.tsv", sep="\t")
a = ax[0, 0]
for i, s in enumerate(sites):
    for j, c in enumerate(classes):
        r = CE[(CE.label == s) & (CE.drug_class == c)]
        if len(r) == 0: continue
        z = r.z_med.iloc[0]; n = r.n_units.iloc[0]
        if not np.isfinite(z): continue
        a.scatter(j, i, s=min(abs(z), 10) * 28 + 6, c="#2a78d6" if z < 0 else "#eb6834", alpha=.85, edgecolors="white", lw=1)
        a.text(j + .28, i, f"{z:.1f}", va="center", fontsize=6.5, color="#444")
a.set_yticks(range(len(sites))); a.set_yticklabels(sites); a.invert_yaxis()
a.set_xticks(range(len(classes))); a.set_xticklabels(classes, rotation=30)
for i in [4.5, 9.5, 11.5]: a.axhline(i, color="#bbb", lw=.8)
a.set_title("A  Median z per drug class (study×agent×cell units)\nblue = decrease, size ∝ |z|; rows: candidates | direct mTORC1 | S6K | AKT", loc="left")
a.grid(False)

a = ax[0, 1]
show = ["LAMTOR1 S98", "TBC1D5 S730", "4EBP1 T46", "4EBP1 S65", "PRAS40 S183", "ULK1 S758", "LARP1 S766", "EIF4B S422", "PRAS40 T246", "GSK3B S9"]
for k, (con, off) in enumerate([("S6Ki/cat_mTORi", -.18), ("AKTi/cat_mTORi", .18)]):
    m = M[(M.contrast == con) & M.label.isin(show) & M.ratio.notna()]
    for i, s in enumerate(show):
        v = m[m.label == s].ratio.clip(-1, 1.5)
        g = m[m.label == s].group
        c = "#eb6834" if k == 0 else "#1baf7a"
        if len(v) == 0: continue
        a.scatter(np.full(len(v), i + off), v, s=26, c=c, edgecolors="white", lw=1, label=con if i == 0 else None)
a.axhline(0, color="#888", lw=.8); a.axhline(1, color="#888", lw=.6, ls="--")
a.set_xticks(range(len(show))); a.set_xticklabels(show, rotation=40, ha="right")
a.set_ylabel("z(inhibitor) / z(catalytic mTORi), same study & cell line")
a.legend(frameon=False, loc="upper left")
a.set_title("B  Matched ratio z(S6Ki or AKTi)/z(mTORi), where mTORi z<−2\n0 = insensitive (direct-like); 1 = as sensitive as to mTORi", loc="left")

a = ax[1, 0]
lab = ["LAMTOR1 S97", "PRAS40 S183", "LARP1 S766", "RRAGC T394", "PRAS40 T246", "TSC2 S939", "GSK3B S9", "RPS6 S244"]
drugs = [("DS0009_AZD80550nM_0", "AZD8055"), ("DS0014_Dactolisib0nM_0", "dactolisib"), ("DS0017_Pictilisib0nM_0", "pictilisib")]
grp = {"LAMTOR1 S97": "cand", "PRAS40 S183": "direct_mTORC1", "LARP1 S766": "direct_mTORC1", "RRAGC T394": "direct_mTORC1",
       "PRAS40 T246": "AKT", "TSC2 S939": "AKT", "GSK3B S9": "AKT", "RPS6 S244": "S6K"}
for j, (e, nm) in enumerate(drugs):
    for i, l in enumerate(lab):
        r = D[(D.exp_id == e) & (D.label == l) & D.sig & (D.delta < 0)]
        if len(r) == 0: continue
        a.scatter(j + (i - 3.5) * .07, r.pEC50.iloc[0], s=40, c=COL[grp[l]], edgecolors="white", lw=1, marker="D" if l == "LAMTOR1 S97" else "o")
        if l == "LAMTOR1 S97": a.annotate("LAMTOR1 S97/98", (j + (i - 3.5) * .07, r.pEC50.iloc[0]), xytext=(4, -10), textcoords="offset points", fontsize=7)
a.set_xticks(range(3)); a.set_xticklabels([d[1] for d in drugs]); a.set_ylabel("pEC50 (decryptM, A549, 1 h)")
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([], [], marker="o", ls="", color=COL[k], label=k) for k in ["direct_mTORC1", "AKT", "S6K"]] +
         [Line2D([], [], marker="D", ls="", color=COL["cand"], label="LAMTOR1 S97 (same peptide as S98)")], frameon=False, fontsize=7, loc="center left")
a.set_title("C  Dose: pictilisib separates AKT-driven (pEC50≈6.8–7) from\nmTOR-direct sites (≈5.4–6.4, needs mTOR off-target)", loc="left")

a = ax[1, 1]
cs = ["LAMTOR1 S98", "TBC1D5 S730", "PRKAA1 S356", "4EBP1 S65", "PRAS40 S183", "ULK1 S758", "LARP1 S766", "RPS6 S240", "EIF4B S422", "GSK3B S9", "TSC2 S939"]
rd = [("4EBP1_T37T46", "#2a78d6"), ("RPS6_S235S236", "#eb6834"), ("PRAS40_T246", "#1baf7a")]
for k, (r, c) in enumerate(rd):
    v = [CP[(CP.site == s) & (CP.readout == r) & (CP.kind == "partial_mTOR")].meta_z.squeeze() if len(CP[(CP.site == s) & (CP.readout == r) & (CP.kind == "partial_mTOR")]) else np.nan for s in cs]
    a.barh(np.arange(len(cs)) + (k - 1) * .27, np.clip(np.array(v, float), -6, 20), height=.25, color=c, label=r)
a.set_yticks(range(len(cs))); a.set_yticklabels(cs); a.invert_yaxis(); a.axvline(0, color="#888", lw=.8)
a.set_xlabel("CPTAC partial co-variation, Stouffer z (protein-adjusted; clipped at 20)")
a.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(.5, -.12), ncol=3); a.set_title("D  CPTAC partial co-variation with three readouts\n(RRAGC S15 not quantified in CPTAC)", loc="left")
fig.tight_layout(); fig.savefig(f"{F}/fig_T1_mtor.png", dpi=170); plt.close(fig)

# ---------------- T2 ----------------
P = pd.read_csv(f"{T}/t2_t3_percentiles.tsv", sep="\t", dtype={"study_pmid": str})
C3 = pd.read_csv(f"{T}/t3_mapk_concordance.tsv", sep="\t", dtype={"study_pmid": str})
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
a = ax[0]
k = U[(U.label == "KPNA4 S60") & U.drug_class.isin(["CK2i", "MEKi", "ERKi", "RAFi", "EGFRi", "cat_mTORi", "AKTi", "S6Ki"])]
order = ["CK2i", "MEKi", "ERKi", "RAFi", "EGFRi", "cat_mTORi", "AKTi", "S6Ki"]
for i, c in enumerate(order):
    v = k[k.drug_class == c]
    cc = ["#2a78d6" if c == "CK2i" else "#9a9a94"] * len(v)
    if len(v) == 0: continue
    a.scatter(np.full(len(v), i) + np.random.default_rng(1).uniform(-.15, .15, len(v)), v.z, s=24, c=cc, edgecolors="white", lw=1)
    for _, r in v[v.z < -2].iterrows(): a.annotate(f"{r.agents} {r.cell_line_norm}", (i, r.z), xytext=(5, 0), textcoords="offset points", fontsize=6)
a.axhline(0, color="#888", lw=.8); a.set_xticks(range(len(order))); a.set_xticklabels(order, rotation=30)
a.set_ylabel("KPNA4 S60 robust z"); a.set_title("A  KPNA4 S60 by inhibitor class (each dot = study×agent×cell)", loc="left")
a = ax[1]
t = C3.dropna(subset=["KPNA4_S60"])
a.scatter(t.ERK_loop.clip(-15, 12), t.KPNA4_S60, s=24, c=np.where(t.pathway_down, "#eda100", "#9a9a94"), edgecolors="white", lw=1)
a.axhline(0, color="#888", lw=.8); a.axvline(0, color="#888", lw=.8)
a.set_xlabel("ERK1/2 pT/pY z (clipped)"); a.set_ylabel("KPNA4 S60 z")
a.set_title("B  MEK/ERK/RAF/EGFR inhibitor units: KPNA4 S60 does not\nfollow ERK (yellow = pathway verified down; 0/25 with z<−2)", loc="left")
fig.tight_layout(); fig.savefig(f"{F}/fig_T2_kpna4.png", dpi=170); plt.close(fig)

# ---------------- T3 ----------------
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
a = ax[0]
t = C3.dropna(subset=["EPN1_S435", "ERK_loop"])
cmap = {"MEKi": "#2a78d6", "ERKi": "#eb6834", "RAFi": "#1baf7a", "EGFRi": "#eda100"}
for c, g in t.groupby("drug_class"):
    a.scatter(g.ERK_loop.clip(-15, 12), g.EPN1_S435, s=26, c=cmap[c], edgecolors="white", lw=1, label=c)
a.axhline(0, color="#888", lw=.8); a.axvline(0, color="#888", lw=.8)
rho = t[["ERK_loop", "EPN1_S435"]].corr("spearman").iloc[0, 1]
a.set_xlabel("ERK1/2 pT/pY z (clipped)"); a.set_ylabel("EPN1 S435 z"); a.legend(frameon=False, fontsize=7)
a.set_title(f"A  Units: EPN1 S435 vs ERK activation-loop\nSpearman ρ = {rho:.2f} (n={len(t)})", loc="left")
a = ax[1]
exps = [("DS0018_Refametinib0nM_0", "refametinib A549"), ("DS0045_Lapatinib0nM_0", "lapatinib MDA-MB-175"), ("DS0030_Lapatinib0nM_0", "lapatinib BT474")]
labs = [("EPN1 S435", "#222222"), ("ERK1 Y204", "#eda100"), ("ERK2 Y187", "#eda100"), ("ERK1 T202", "#eda100"), ("RSK1 S380", "#eb6834")]
for j, (e, nm) in enumerate(exps):
    for i, (l, c) in enumerate(labs):
        r = D[(D.exp_id == e) & (D.label == l) & (D.delta < 0) & (D.R2 > .8)]
        if len(r) == 0: continue
        a.scatter(j + (i - 2) * .08, r.pEC50.iloc[0], s=40, c=c, marker="D" if l.startswith("EPN1") else "o", edgecolors="white", lw=1,
                  alpha=1 if r.sig.iloc[0] else .5)
a.set_xticks(range(3)); a.set_xticklabels([e[1] for e in exps], rotation=15); a.set_ylabel("pEC50 (R²>0.8; faded = q>0.05)")
a.legend(handles=[Line2D([], [], marker="D", ls="", color="#222", label="EPN1 S435"), Line2D([], [], marker="o", ls="", color="#eda100", label="ERK1/2 loop"),
                  Line2D([], [], marker="o", ls="", color="#eb6834", label="RSK1 S380")], frameon=False, fontsize=7)
a.set_title("B  decryptM potency: EPN1 S435 tracks ERK", loc="left")
a = ax[2]
cs = ["EPN1 S435", "RSK1 S380", "KPNA4 S60", "GSK3B S9"]
rd = [("ERK_pTpY", "#eda100"), ("CK2_comp", "#2a78d6"), ("CDK_comp", "#9a9a94")]
for k2, (r, c) in enumerate(rd):
    v = [(lambda q: q.meta_z.iloc[0] if len(q) else np.nan)(CP[(CP.site == s) & (CP.readout == r) & (CP.kind == "partial_MAPK")]) for s in cs]
    a.barh(np.arange(len(cs)) + (k2 - 1) * .27, np.array(v, float), height=.25, color=c, label=r)
a.set_yticks(range(len(cs))); a.set_yticklabels(cs); a.invert_yaxis(); a.axvline(0, color="#888", lw=.8)
a.legend(frameon=False, fontsize=7); a.set_xlabel("CPTAC partial Stouffer z")
a.set_title("C  Tumours: partial co-variation (ERK | CK2 | CDK)", loc="left")
fig.tight_layout(); fig.savefig(f"{F}/fig_T3_epn1.png", dpi=170); plt.close(fig)
print("ok")
