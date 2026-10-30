"""Audit H figure: (A) H1 split-half by split scheme; (B) H3 coupling under adjustments; (C) H4 per-study LOLO deltas."""
import os; os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
O = "/workspace/results/audit_H"; B, OR, AQ, INK, MUT = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUT,
                     "axes.labelcolor": INK, "xtick.color": MUT, "ytick.color": MUT})
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1.2]))
H = pd.read_csv(f"{O}/h1_lab_splits.tsv", sep="\t"); H = H[H.group == "ALL"]
N = pd.read_csv(f"{O}/h1_noise_matched_random.tsv", sep="\t"); N = N[(N.group == "ALL") & N.level.isin(["study", "lab"])]
sch = [("study_split (I28 replicate)", "study split\n(I28)"), ("lab_split", "lab split"), ("Cutillas_vs_rest (fixed)", "Cutillas vs\nrest"),
       ("study_split_drop_Cutillas", "study split\n−Cutillas"), ("study_split_drop_Ahn", "study split\n−Ahn")]
x = np.arange(len(sch)); w = .27
for k, (s, c, lab) in enumerate([("substrate", B, "substrates"), ("random", OR, "matched random (I28)")]):
    v = [H[(H.analysis == a) & (H.set == s)].rho.iloc[0] for a, _ in sch]
    ax[0].bar(x + (k - .5) * w, v, w * .92, color=c, label=lab)
nm = N[N.set == "random_noise_matched"].set_index("level").rho
ax[0].bar([0 + 1.5 * w, 1 + 1.5 * w], [nm["study"], nm["lab"]], w * .92, color=AQ, label="random, + noise/dynamism/coverage matched")
ax[0].set_xticks(x + .1); ax[0].set_xticklabels([l for _, l in sch]); ax[0].set_ylabel("split-half Spearman (within family)")
ax[0].set_title("A  H1: reproducibility shrinks across labs but persists", loc="left", color=INK); ax[0].legend(frameon=False, fontsize=8, loc="upper left"); ax[0].set_ylim(0, .68)
S = pd.read_csv(f"{O}/h3_summary.tsv", sep="\t"); S = S[(S.subset == "all") & (S.outcome == "coupling")]
x = np.arange(3)
ax[1].bar(x - .2, S.site_median_fast, .38, color=B, label="fast substrates"); ax[1].bar(x + .2, S.site_median_stable, .38, color=OR, label="stable substrates")
for i, r in enumerate(S.itertuples()):
    ax[1].text(i, max(r.site_median_fast, r.site_median_stable) + .01, f"adj Δ {r.adj_diff:.3f}\n[{r.adj_lo:.3f},{r.adj_hi:.3f}]", ha="center", fontsize=7.5, color=MUT)
ax[1].set_xticks(x); ax[1].set_xticklabels(["I28\n(protein-resid.)", "+ ischaemia", "+ ischaemia\n+ 5 PCs"])
ax[1].set_ylabel("median site coupling to family footprint (ρ)"); ax[1].set_ylim(0, .32)
ax[1].set_title("B  H3: fast–stable gap roughly halves", loc="left", color=INK); ax[1].legend(frameon=False, fontsize=8, loc="upper right")
R = pd.read_csv(f"{O}/h4_activity_conditions_lolo.tsv", sep="\t", dtype={"study": str})
P = R.groupby(["study", "condition", "vbase", "target_edited"]).auroc.mean().unstack("vbase").reset_index()
P = P[P.target_edited].dropna(subset=["base"]); P["d"] = P.keep_fast_only - P.base; P["dr"] = P.rand_keepn - P.base
s = P.groupby("study").agg(n=("d", "size"), d=("d", "mean"), dr=("dr", "mean")).sort_values("d")
y = np.arange(len(s))
ax[2].hlines(y, s.dr, s.d, color="#c8c7c2", lw=1.5)
ax[2].scatter(s.dr, y, s=10 + 2 * s.n, color=OR, label="size-matched random drop", zorder=3)
ax[2].scatter(s.d, y, s=10 + 2 * s.n, color=B, label="keep fast edges only", zorder=3)
ax[2].axvline(0, color=MUT, lw=.8); ax[2].set_yticks(y); ax[2].set_yticklabels([f"{i} (n={n})" for i, n in zip(s.index, s.n)], fontsize=7)
ax[2].set_xlabel("ΔAUROC vs base (leave-one-LAB-out turnover)")
ax[2].set_title("C  H4: gain from few studies; decryptM ≈ 0", loc="left", color=INK); ax[2].legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1, .12), markerscale=.5)
fig.tight_layout(); fig.savefig(f"{O}/fig_audit_H.png", dpi=170)
