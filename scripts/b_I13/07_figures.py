import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = "/workspace/results/b_I13"; FG = f"{R}/figures"
BLUE, ORANGE, AQUA, GREY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#c3c2b7", "#0b0b0b", "#898781"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": GREY, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
# Fig 1: perturbation univariate effects heatmap
u = pd.read_csv(f"{R}/tables/t1_univariate_logit.tsv", sep="\t").pivot_table(index="feature", columns="cond", values="beta")
fig, ax = plt.subplots(figsize=(5.5, 6)); v = 1.3
im = ax.imshow(u.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
ax.set_xticks(range(u.shape[1]), u.columns, rotation=40, ha="right"); ax.set_yticks(range(u.shape[0]), u.index)
for i in range(u.shape[0]):
    for j in range(u.shape[1]): ax.text(j, i, f"{u.values[i,j]:.2f}", ha="center", va="center", fontsize=6.5, color=INK)
fig.colorbar(im, ax=ax, shrink=0.6, label="log-odds per SD (site up, z>=2)")
ax.set_title("Compendium phosphatase inhibitors: feature → site increase\n(OA, CA = PP2A inhib., MCF7, 1 study; PV = pervanadate)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{FG}/fig1_perturbation_features.png", dpi=150); plt.close()
# Fig 2: per-cohort forest for PPP2R2A exposures (pc5)
s = pd.read_csv(f"{R}/tables/t2_stage2_per_cohort_pc5.tsv", sep="\t")
tests = ["bipartite_basic|ST", "B55_lit|ST", "isT|ST", "P1|ST"]
fig, axes = plt.subplots(1, 4, figsize=(11, 3.6), sharey=True)
for ax, t in zip(axes, tests):
    for k, (ex, col, off) in enumerate([("PPP2R2A_cnvloss", BLUE, -0.15), ("PPP2R2A_protlow", ORANGE, 0.15)]):
        d = s[(s.exposure == ex) & (s.test == t)].sort_values("cohort")
        y = np.arange(len(d)) + off
        ax.scatter(d.z_emp, y, color=col, s=22, label=ex if t == tests[0] else None, zorder=3)
    ax.axvline(0, color=GREY, lw=1); ax.axvspan(-1.96, 1.96, color="#e1e0d9", alpha=0.4, lw=0)
    ax.set_title(t.replace("|ST", " (pS/pT)"), fontsize=9); ax.set_xlabel("empirical z vs random-gene/protein null")
axes[0].set_yticks(range(10), sorted(s.cohort.unique())); axes[0].legend(frameon=False, fontsize=7, loc="lower left")
fig.suptitle("PPP2R2A loss (CNV) and low PPP2R2A protein vs B55-rule site features, per CPTAC cohort (protein-adjusted sites; purity, E2F/G2M, CNV burden, immune/stroma, 5 proteome PCs)", fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{FG}/fig2_PPP2R2A_forest.png", dpi=150); plt.close()
# Fig 3: rank vs 1000 random proteins
S2 = pd.read_parquet(f"{R}/cache/stage2_all_pc5.parquet")
d = S2[S2.test == "bipartite_basic|ST"]; npr = d[d.exposure.str.startswith("null_prot")]
st = npr.groupby("cohort").coef.agg(["mean", "std"]); d = d.join(st, on="cohort"); d["z"] = (d.coef - d["mean"]) / d["std"]
M = d.groupby("exposure").z.agg(lambda z: z.sum() / np.sqrt(len(z)))
nul = M[M.index.str.startswith("null_prot")]
fig, ax = plt.subplots(figsize=(5.5, 3.2))
ax.hist(nul, bins=40, color=GREY, edgecolor="#fcfcfb")
for ex, col in [("PPP2R2A_protlow", ORANGE), ("PPP2R2A_cnvloss", BLUE), ("PPP1CB_protlow", AQUA)]:
    if ex in M: ax.axvline(M[ex], color=col, lw=2); ax.text(M[ex], ax.get_ylim()[1] * 0.9, ex, rotation=90, ha="right", va="top", fontsize=7, color=INK)
ax.set_xlabel("Stouffer z across 10 cohorts (bipartite-basic pS/pT sites up when exposure low)"); ax.set_ylabel("random proteins (n=1000)")
ax.set_title("Strict null: same random protein used in every cohort", fontsize=9)
fig.tight_layout(); fig.savefig(f"{FG}/fig3_protein_null_rank.png", dpi=150); plt.close()
# Fig 4: kinase activity vs PPP2R2A CNV loss: vs other-chromosome null and vs 8p null
k = pd.read_csv(f"{R}/tables/t3_kinase_activity_meta_pc5.tsv", sep="\t")
k = k[(k.variant == "tc") & (k.exposure == "PPP2R2A_cnvloss")].sort_values("stouffer")
fig, ax = plt.subplots(figsize=(5.5, 5)); y = np.arange(len(k))
ax.scatter(k.stouffer, y, color=BLUE, s=22, label="vs random genes, other chromosomes", zorder=3)
ax.scatter(k.stouffer_vs8p, y, color=ORANGE, s=22, marker="D", label="vs random 8p genes (same arm)", zorder=3)
ax.set_yticks(y, k.kinase); ax.axvline(0, color=GREY); ax.axvspan(-1.96, 1.96, color="#e1e0d9", alpha=0.4, lw=0)
ax.set_xlabel("Stouffer z (core-v3 kinase activity ~ PPP2R2A CNV loss)"); ax.legend(frameon=False, fontsize=7, loc="lower right")
ax.set_title("Kinase activity vs PPP2R2A CNV loss (10 cohorts)\nHIPK2/DYRK1A dip vanishes vs same-arm null", fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{FG}/fig4_kinase_activity_PPP2R2A.png", dpi=150); plt.close()
print("ok")
