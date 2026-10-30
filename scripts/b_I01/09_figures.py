"""Figures for b_I01."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
OUT = Path("/workspace/results/b_I01"); TB = OUT / "tables"; FG = OUT / "figures"
BLUE, ORANGE, AQUA, YELLOW, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8984"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#52514e",
                     "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e"})

# Fig 1: family echo and specific-signal coherence
T = pd.read_csv(TB / "specificity_per_kinase.tsv", sep="\t", index_col=0)
T10 = pd.read_csv(TB / "specificity_per_kinase_pc10.tsv", sep="\t", index_col=0)
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
for j, (col, lab) in enumerate([("r2_full_conf", "all predicted substrates"), ("r2_spec_conf", "kinase-specific substrates")]):
    for k, (dk, c) in enumerate([(True, ORANGE), (False, BLUE)]):
        v = T.loc[T.dark == dk, col].dropna()
        ax[0].boxplot(v, positions=[j * 2.6 + k], widths=0.7, patch_artist=True, showfliers=False,
                      boxprops=dict(facecolor=c, alpha=.35, edgecolor=c), medianprops=dict(color=c, lw=2), whiskerprops=dict(color=c), capprops=dict(color=c))
ax[0].set_xticks([0.5, 3.1]); ax[0].set_xticklabels(["all predicted\nsubstrates", "kinase-specific\nsubstrates"])
ax[0].set_ylabel("R² with family + motif-neighbour score"); ax[0].set_title("Family echo (orange = dark, blue = studied)", fontsize=9)
for a, D, t in [(ax[1], T, "Coherence beyond family"), (ax[2], T10, "…and beyond 10 global kinome PCs")]:
    d = D.dropna(subset=["rel_obs"])
    for dk, c in [(False, BLUE), (True, ORANGE)]:
        x = d[d.dark == dk]
        a.scatter(x.rel_null_mean, x.rel_obs, s=10, color=c, alpha=.6, label="dark" if dk else "studied", edgecolor="white", lw=.3)
    lim = [min(d.rel_null_mean.min(), d.rel_obs.min()) - .02, max(d.rel_obs.max(), .1) + .02]
    a.plot(lim, lim, color=GREY, lw=1, ls="--"); a.set_xlabel("split-half r, matched random site sets")
    a.set_ylabel("split-half r, specific sites (protein-split)"); a.set_title(t, fontsize=9)
for n in ["VRK2", "TLK1", "MAPK15", "DSTYK", "CLK4", "DYRK3"]:
    if n in T10.index and pd.notna(T10.loc[n, "rel_obs"]):
        ax[2].annotate(n, (T10.loc[n, "rel_null_mean"], T10.loc[n, "rel_obs"]), fontsize=7, color="#0b0b0b", xytext=(3, 2), textcoords="offset points")
ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FG / "fig1_family_echo_coherence.png", dpi=180); plt.close(fig)

# Fig 2: inhibitor positive control (conditions where drug binds k but none of k's family/neighbours)
fig, ax = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True)
for a, tag, t in [(ax[0], "", "family-partialled residual"), (ax[1], "_pc10", "+ 10 kinome PCs removed")]:
    pc = pd.read_csv(TB / f"inhibitor_control{tag}.tsv", sep="\t")
    pc = pc[(pc.target_def == "profiled") & pc.conf_hit_excluded]
    for i, (dk, lab) in enumerate([(False, "studied"), (True, "dark")]):
        g = pc[pc.dark == dk]
        for j, (m, c) in enumerate([("full", GREY), ("conf", BLUE), ("resid", ORANGE)]):
            v = g["auc_" + m].dropna()
            pos = i * 4 + j
            a.boxplot(v, positions=[pos], widths=.7, patch_artist=True, showfliers=False, boxprops=dict(facecolor=c, alpha=.35, edgecolor=c),
                      medianprops=dict(color=c, lw=2), whiskerprops=dict(color=c), capprops=dict(color=c))
            a.text(pos, 1.02, f"n={len(v)}", ha="center", fontsize=6.5, color="#52514e")
    a.axhline(.5, color=GREY, lw=1, ls="--"); a.set_xticks([1, 5]); a.set_xticklabels(["studied kinases", "dark kinases"])
    a.set_title(t, fontsize=9); a.set_ylim(0, 1.08)
ax[0].set_ylabel("AUROC: score lower when drug binds kinase\n(grey all-substrate, blue family, orange specific)")
fig.tight_layout(); fig.savefig(FG / "fig2_inhibitor_control.png", dpi=180); plt.close(fig)

# Fig 3: orthogonal validation
R = pd.read_csv(TB / "validation_auroc.tsv", sep="\t")
labs = ["depmap_noprox", "ppi_exp", "ppi_lit", "gtex", "hpa", "compartments"]
meths = [("full", GREY, "all substrates"), ("conf", BLUE, "family/neighbours"), ("resid", ORANGE, "specific residual"), ("family_label", AQUA, "family-label baseline")]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
for a, qs, t in [(ax[0], "dark_reliable", "26 dark kinases with coherent specific signal"), (ax[1], "nondark", "308 studied kinases")]:
    d = R[R.queries == qs]
    for i, l in enumerate(labs):
        for j, (m, c, _) in enumerate(meths):
            r = d[(d.label == l) & (d.method == m)]
            if len(r) == 0: continue
            x = i + (j - 1.5) * .19
            a.bar(x, r.auroc.iloc[0] - .5, bottom=.5, width=.17, color=c)
            if m != "family_label" and r.p_perm.iloc[0] <= 0.05 and r.auroc.iloc[0] > .5:
                a.text(x, max(r.auroc.iloc[0], .5) + .005, "*", ha="center", fontsize=9)
    a.axhline(.5, color=GREY, lw=1); a.set_xticks(range(len(labs)))
    a.set_xticklabels(["DepMap\nco-dep", "PPI\nscreens", "PPI\nliterature", "GTEx\nco-expr", "HPA\nco-loc", "COMPART-\nMENTS"], fontsize=7.5)
    a.set_title(t, fontsize=9); a.set_ylim(.35, .85)
ax[0].set_ylabel("AUROC (non-family pairs)")
ax[1].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in meths], labels=[l for *_, l in meths], frameon=False, fontsize=7.5, loc="upper right")
fig.text(.01, .01, "* row-permutation p ≤ 0.05 (300 perms), AUROC > 0.5 only; DepMap excludes same-cytoband pairs", fontsize=7, color="#52514e")
fig.tight_layout(); fig.savefig(FG / "fig3_orthogonal_validation.png", dpi=180); plt.close(fig)
print("ok")
