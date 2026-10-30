"""I28 figures. fig1 atlas+reliability; fig2 predictors; fig3 consequences (CPTAC, activity)."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); from lib28 import *
BLUE, ORANGE, GRAY, INK, MUTED = "#2a78d6", "#eb6834", "#9a9a94", "#262624", "#6b6a64"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6})
T = f"{R28}/tables"
A = pd.read_csv(f"{T}/turnover_atlas.tsv", sep="\t"); L = pd.read_parquet(f"{C28}/L_units.parquet")
AR = classify_units(L[L.set == "random"])
GS = pd.read_csv(f"{T}/group_summary.tsv", sep="\t"); RL = pd.read_csv(f"{T}/reliability.tsv", sep="\t")
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for i, g in enumerate(MAIN):
    a = A[(A.group == g) & (A.n_units >= 2)].R.clip(-8, 4); r = AR[(AR.group == g) & (AR.n_units >= 2)].R.clip(-8, 4)
    for x, v, c in [(i - 0.18, r, GRAY), (i + 0.18, a, BLUE)]:
        bp = ax[0].boxplot(v, positions=[x], widths=0.3, patch_artist=True, showfliers=False)
        for p in bp["boxes"]: p.set(facecolor=c, alpha=0.6, edgecolor=c)
        for k in ["medians", "whiskers", "caps"]:
            for p in bp[k]: p.set(color=c)
ax[0].axhline(-1, ls="--", color=MUTED, lw=0.8)
ax[0].set_xticks(range(len(MAIN))); ax[0].set_xticklabels([f"{g}\n(n={len(A[(A.group == g) & (A.n_units >= 2)])})" for g in MAIN], fontsize=7)
ax[0].set_ylabel("responsiveness R (mean robust z, verified ≤2 h inhibition)"); ax[0].set_title("A  Substrates (blue) vs matched random (gray)", loc="left")
gs = GS[GS.group.isin(MAIN)].set_index("group").loc[MAIN]
x = np.arange(len(MAIN))
ax[1].bar(x - 0.2, gs.frac_fast_random, 0.38, color=GRAY, label="matched random")
ax[1].bar(x + 0.2, gs.frac_fast, 0.38, color=BLUE, label="confident substrates")
ax[1].errorbar(x + 0.2, gs.frac_fast, yerr=[gs.frac_fast - gs.frac_fast_lo, gs.frac_fast_hi - gs.frac_fast], fmt="none", ecolor=INK, lw=1)
ax[1].set_xticks(x); ax[1].set_xticklabels(MAIN, fontsize=7); ax[1].set_ylabel("fraction fast"); ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title("B  Fraction 'fast' (95% CI)", loc="left")
rl = RL[(RL.group.isin(["ALL"] + MAIN)) & (RL.level == "study") & RL.set.isin(["substrate", "random"])]
gg = ["ALL"] + [g for g in MAIN if g in set(rl.group)]
for j, (s, c) in enumerate([("random", GRAY), ("substrate", BLUE)]):
    d = rl[rl.set == s].set_index("group").reindex(gg)
    ax[2].bar(np.arange(len(gg)) + (j - 0.5) * 0.38, d.rho, 0.38, color=c, label=s)
    ax[2].errorbar(np.arange(len(gg)) + (j - 0.5) * 0.38, d.rho, yerr=[d.rho - d.rho_lo, d.rho_hi - d.rho], fmt="none", ecolor=INK, lw=1)
ax[2].set_xticks(range(len(gg))); ax[2].set_xticklabels(gg, fontsize=7); ax[2].set_ylabel("split-half Spearman (across studies)")
ax[2].set_title("C  Cross-study reproducibility", loc="left"); ax[2].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{R28}/figures/fig1_atlas_reliability.png", dpi=160); plt.close()

U = pd.read_csv(f"{T}/predictors_univariate.tsv", sep="\t").sort_values("OR_per_SD")
U = U[U.feature != "lability_other"]
LS = pd.read_csv(f"{T}/predictors_loso_summary.tsv", sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(12, 7), gridspec_kw={"width_ratios": [1.4, 1]})
y = np.arange(len(U)); c = np.where(U.q < 0.05, BLUE, GRAY)
ax[0].hlines(y, U.lo, U.hi, color=c, lw=1.5); ax[0].scatter(U.OR_per_SD, y, color=c, s=18, zorder=3)
ax[0].axvline(1, color=MUTED, lw=0.8); ax[0].set_yticks(y); ax[0].set_yticklabels(U.feature, fontsize=7); ax[0].set_xscale("log")
ax[0].set_xlabel("OR per SD for fast vs stable (family FE + coverage; protein-clustered 95% CI)")
ax[0].set_title("A  Single-feature associations (blue: FDR<5%)", loc="left")
order = ["pred_motif_only", "pred_features", "pred_features+lability", "R_other"]
lab = {"pred_motif_only": "motif only", "pred_features": "intrinsic features", "pred_features+lability": "features + other-inhibitor lability", "R_other": "same site, other studies"}
ls = LS.set_index("predictor").loc[order]
ax[1].barh(range(4), ls.auroc_wmean - 0.5, left=0.5, color=[GRAY, BLUE, BLUE, ORANGE])
for i, v in enumerate(ls.auroc_wmean): ax[1].text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=8, color=INK)
ax[1].set_yticks(range(4)); ax[1].set_yticklabels([lab[o] for o in order], fontsize=8); ax[1].set_xlim(0.45, 0.8)
ax[1].set_xlabel("AUROC for falling (z≤-1) in held-out study\n(sites also held out for feature models)")
ax[1].set_title("B  Leave-one-study-out prediction", loc="left")
fig.tight_layout(); fig.savefig(f"{R28}/figures/fig2_predictors.png", dpi=160); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
if os.path.exists(f"{T}/cptac_site_stats.tsv"):
    D = pd.read_csv(f"{T}/cptac_site_stats.tsv", sep="\t")
    d = D.groupby(["kind", "cls", "group", "site"]).coupling.mean().reset_index()
    cats = [("random", None, "random", GRAY), ("substrate", "stable", "stable", "#b7d3f6"), ("substrate", "intermediate", "intermediate", "#6da7ec"), ("substrate", "fast", "fast", BLUE)]
    for i, (k, cl_, lb, col) in enumerate(cats):
        v = d[(d.kind == k) & ((d.cls == cl_) if cl_ else True)].coupling.dropna()
        bp = ax[0].boxplot(v, positions=[i], widths=0.5, patch_artist=True, showfliers=False)
        for p in bp["boxes"]: p.set(facecolor=col, edgecolor=MUTED)
        ax[0].text(i, ax[0].get_ylim()[1] if False else 0.75, f"n={len(v)}", ha="center", fontsize=7, color=MUTED)
    ax[0].set_xticks(range(4)); ax[0].set_xticklabels([c[2] for c in cats]); ax[0].set_ylabel("Spearman with LOO family footprint\n(protein-residualised, mean over cohorts)")
    ax[0].set_title("A  CPTAC tumours: coupling to family activity", loc="left")
if os.path.exists(f"{T}/activity_turnover_summary.tsv"):
    O = pd.read_csv(f"{T}/activity_turnover_summary.tsv", sep="\t")
    o = O[O.subset == "target_edited"].set_index("variant")
    vs = ["drop_stable", "rand_drop", "keep_fast_only", "weighted", "rand_weight"]; o = o.reindex(vs)
    cols = [BLUE, GRAY, BLUE, BLUE, GRAY]
    ax[1].barh(range(len(vs)), o.mean_delta_study, color=cols); ax[1].errorbar(o.mean_delta_study, range(len(vs)), xerr=[o.mean_delta_study - o.lo, o.hi - o.mean_delta_study], fmt="none", ecolor=INK)
    ax[1].axvline(0, color=MUTED, lw=0.8); ax[1].set_yticks(range(len(vs))); ax[1].set_yticklabels(vs)
    ax[1].set_xlabel("ΔAUROC vs base PSP scorer (study-mean, 95% bootstrap CI)\nturnover from other studies (LOSO)")
    ax[1].set_title("B  Turnover-aware library, held-out inhibitor studies", loc="left")
fig.tight_layout(); fig.savefig(f"{R28}/figures/fig3_consequences.png", dpi=160); plt.close()
print("figs done")
