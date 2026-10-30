"""fig4: Finding K2 - perturbation-learned programmes in human tissue. Re-plots saved tables only."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from style import *

FP = [("FP:KI_MEK", "MEK inhibitor"), ("FP:EGF", "EGF"), ("FP:TNF", "TNF"), ("FP:DNA_damage", "DNA damage")]
t1 = pd.read_csv(f"{R}/b_Q8/tables/T1_coherence.tsv", sep="\t")
t1 = t1[(t1.cond == "padj_pc3") & (t1.design == "paired")]
t2 = pd.read_csv(f"{R}/b_Q8/tables/T2_coherence_summary.tsv", sep="\t")
t2 = t2[t2.cond == "padj_pc3"].set_index("set")
t7 = pd.read_csv(f"{R}/b_Q8/tables/T7_edge_global_preservation.tsv", sep="\t")
t7 = t7[(t7.cond == "padj") & (t7.label == "edge")].copy()
vb = pd.read_csv(f"{R}/val_B/T3_coherence.tsv", sep="\t")

fig = plt.figure(figsize=(7.2, 5.8))
axA = fig.add_axes([0.10, 0.565, 0.50, 0.36])
axB = fig.add_axes([0.735, 0.565, 0.235, 0.36])
axC = fig.add_axes([0.10, 0.075, 0.885, 0.33])

# ---------- (a) coherence NAT vs tumour ----------
ax = axA
ygrid(ax)
dx = 0.2
cohorts = sorted(t1.cohort.unique())
for i, (s, name) in enumerate(FP):
    d = t1[t1.set == s].set_index("cohort").loc[cohorts]
    col = PINK if s == "FP:DNA_damage" else GREEN
    for c, r in d.iterrows():
        ax.plot([i - dx, i + dx], [r.coh_N, r.coh_T], color=col, lw=0.8, alpha=0.6, zorder=2)
    ax.scatter(np.full(len(d), i - dx), d.coh_N, s=16, color=col, edgecolor="white", lw=0.5, zorder=3)
    ax.scatter(np.full(len(d), i + dx), d.coh_T, s=16, color=col, edgecolor="white", lw=0.5, zorder=3)
    # matched random-site reference
    ax.plot([i - dx, i + dx], [d.null_N.median(), d.null_T.median()], color=GREY, lw=1.2, zorder=2)
    ax.scatter([i - dx, i + dx], [d.null_N.median(), d.null_T.median()], marker="D", s=18, color=GREY,
               edgecolor="white", lw=0.5, zorder=3)
    zN, zT = t2.loc[s, "z_N_med"], t2.loc[s, "z_T_med"]
    ax.text(i, 0.178, f"z {zN:.0f} → {zT:.0f}", ha="center", fontsize=7.3, color=INK)
ax.text(3, 0.161, f"gain in tumour:\nmeta z {t2.loc['FP:DNA_damage', 'meta_z_perm']:.1f}", ha="center", va="top",
        fontsize=7, color=PINK, linespacing=1.0)
ax.set_xticks(range(4))
ax.set_xticklabels([n for _, n in FP])
ax.tick_params(axis="x", which="major", length=0, pad=13)
ax.set_xticks([i + s_ for i in range(4) for s_ in (-dx, dx)], minor=True)
ax.set_xticklabels(["normal", "tumour"] * 4, minor=True, fontsize=6.6, color=MUTED)
ax.tick_params(axis="x", which="minor", length=2, pad=2)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.005, 0.19)
ax.set_ylabel("coherence\n(mean site-to-site correlation)")
ax.text(1.0, 0.011, "matched random sites (median)", ha="center", fontsize=7, color=MUTED)
ax.set_yticks([0, 0.05, 0.10, 0.15])
ax.set_yticklabels(["0", "0.05", "0.10", "0.15"])

# ---------- (b) sign agreement ----------
ax = axB
xgrid(ax)
order = ["CCRCC", "COAD", "HNSCC", "LSCC", "LUAD", "PDAC", "OV", "UCEC"]
t7 = t7.set_index("cohort").loc[order]
y = np.arange(len(order))[::-1]
for yy, (c, r) in zip(y, t7.iterrows()):
    ax.plot([r.signconc_N, r.signconc_T], [yy, yy], color="#bbbbbb", lw=1.2, zorder=1)
ax.scatter(t7.signconc_N, y, s=26, facecolor="white", edgecolor=GREEN, lw=1.3, zorder=3, label="normal tissue")
ax.scatter(t7.signconc_T, y, s=26, color=GREEN, edgecolor="white", lw=0.5, zorder=3, label="tumour")
main6 = t7.iloc[:6]
mN, mT = main6.signconc_N.median(), main6.signconc_T.median()
ax.axvline(0.5, color=GREY, ls=(0, (3, 2)), lw=0.9)
ax.text(0.515, 0.0, "chance", fontsize=7, color=MUTED, va="center")
ax.set_yticks(y)
ax.set_yticklabels([f"{c}" + (" *" if c in ("OV", "UCEC") else "") for c in order], fontsize=7.8)
ax.set_xlim(0.45, 1.0)
ax.set_ylim(-2.0, len(order) - 0.4)
ax.set_xlabel("fraction of 529 relationships\nwith the cell-line sign")
ax.set_xticks([0.5, 0.75, 1.0])
ax.legend(loc="lower right", bbox_to_anchor=(1.02, -0.01), ncol=2, fontsize=7, handletextpad=0.2,
          columnspacing=0.6, borderaxespad=0.1)
ax.text(0.53, 1.6, f"median of 6\nnormal {mN:.2f}\ntumour {mT:.2f}", fontsize=7, color=INK, va="center",
        linespacing=1.05)

# ---------- (c) ischaemia time ----------
ax = axC
ygrid(ax)
ADJ = [("none", "no adjustment", GREEN), ("measured t + log t", "adjusted for measured\nischaemia time", PURPLE),
       ("footprint ISCH score", "adjusted for proxy\nischaemia score", ORANGE)]
SETS = [("MEK", "MEK inhibitor"), ("EGF", "EGF"), ("TNF", "TNF"), ("DNAdmg", "DNA damage")]
w = 0.29
rng = np.random.default_rng(2)
for i, (s, name) in enumerate(SETS):
    for j, (a, lab, col) in enumerate(ADJ):
        d = vb[(vb.set == s) & (vb.adj == a)]
        assert len(d) == 7, (s, a, len(d))
        x = i + (j - 1) * w
        ax.scatter(x + rng.uniform(-0.05, 0.05, len(d)), d.coherence, s=12, color=col, alpha=0.75,
                   edgecolor="none", zorder=3)
        m = d.coherence.median()
        ax.plot([x - 0.09, x + 0.09], [m, m], color=col, lw=2.2, zorder=4, solid_capstyle="butt")
        ax.text(x + 0.1, m, f".{round(m * 1000):03d}", ha="left", va="center", fontsize=6.8, color=INK, zorder=5,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.3, alpha=0.85))
ax.set_xticks(range(4))
ax.set_xticklabels([n for _, n in SETS])
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(0, 0.27)
ax.set_ylabel("coherence in tumours\n(mean site-to-site correlation)")
h = [Line2D([], [], color=c, marker="o", lw=2.2, ms=4, label=l.replace("\n", " ")) for _, l, c in ADJ]
ax.legend(handles=h, loc="upper right", fontsize=7, ncol=1, handlelength=1.6)
ax.text(0.995, 0.66, "dots: 7 cohorts with recorded times (n = 723 tumours); bar: median",
        transform=ax.transAxes, ha="right", fontsize=6.8, color=MUTED)

ptitle(fig, axA, "a", "Footprints are coherent in normal tissue and tumours", dx=0.075)
ptitle(fig, axB, "b", "Relationships keep their sign", dx=0.095)
ptitle(fig, axC, "c", "Recorded ischaemia time does not explain the MAPK axis", dx=0.075)
save(fig, "fig4_tissue_wiring")
