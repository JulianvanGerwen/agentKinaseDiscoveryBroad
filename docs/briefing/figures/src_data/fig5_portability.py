"""fig5: Finding K3 - primary vs reverse responses. Re-plots saved tables only."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from style import *

m = pd.read_csv(f"{R}/audit_I19/T_fig_mirror.tsv", sep="\t")
lo = pd.read_csv(f"{R}/audit_I10/A3_loso_own.tsv", sep="\t").dropna(subset=["down", "up"])

fig = plt.figure(figsize=(7.2, 5.9))
W, H = 0.37, 0.30
aI = fig.add_axes([0.105, 0.585, W, H])
aS = fig.add_axes([0.605, 0.585, W, H], sharey=aI)
bI = fig.add_axes([0.105, 0.085, W, 0.33])
bS = fig.add_axes([0.605, 0.085, W, 0.33], sharey=bI)

BINS = ["0.5-1", "1-2", "2-99"]
BLAB = ["0.5–1", "1–2", "> 2"]
for ax, panel, prim, title, nst in [(aI, "Kinase inhibitors", "down", "Kinase inhibitors", "41 studies"),
                                    (aS, "Activating stimulations", "up", "Stimulations (growth factors, TNF, ...)", "52 studies")]:
    ygrid(ax)
    ax.axhline(0, color="#999999", lw=0.8)
    for d, off in [("down", -0.06), ("up", 0.06)]:
        s = m[(m.panel == panel) & (m.dir == d)].set_index("bin").loc[BINS]
        isprim = d == prim
        col = GREEN if isprim else ORANGE
        x = np.arange(3) + off
        ax.errorbar(x, s.lor, yerr=[s.lor - s.lo, s.hi - s.lor], fmt="none", ecolor=col, elinewidth=1.1,
                    capsize=0, alpha=0.8, zorder=2)
        ax.plot(x, s.lor, "-o", color=col, lw=2, ms=5.5, mec="white", mew=0.8, zorder=3)
        word = "decreases" if d == "down" else "increases"
        role = "primary" if isprim else "reverse"
        ax.text(x[-1] + 0.17, s.lor.iloc[-1], f"{word}\n({role})", fontsize=7.5, color=INK, va="center",
                ha="left", linespacing=1.0)
        for xi, v in zip(x, s.lor):
            pass
    ax.set_xticks(range(3))
    ax.set_xticklabels(BLAB)
    ax.set_xlim(-0.35, 3.05)
    ax.set_xlabel("size of change in first study (|log2 fold change|)")
    ax.text(0.0, 1.02, f"{title}  ·  {nst}", transform=ax.transAxes, fontsize=8, color=MUTED, ha="left",
            va="bottom")
aI.set_ylabel("replication in another study\n(log odds ratio vs base rate)")
aI.set_ylim(-0.1, 2.35)
plt.setp(aS.get_yticklabels(), visible=False)

# ---------- (b) footprint transfer, down-site vs up-site ----------
for ax, kind, prim, title in [(bI, "inhibitor", "down", "Inhibitors"),
                              (bS, "stimulus/stress", "up", "Stimuli")]:
    d = lo[lo.kind == kind]
    ax.plot([0.15, 1], [0.15, 1], color=GREY, lw=0.9, ls=(0, (3, 2)), zorder=1)
    ax.axhline(0.5, color=GRID, lw=0.8, zorder=0)
    ax.axvline(0.5, color=GRID, lw=0.8, zorder=0)
    pw = (d.down > d.up) if prim == "down" else (d.up > d.down)
    ax.scatter(d.down[pw], d.up[pw], s=30, color=GREEN, edgecolor="white", lw=0.6, zorder=3)
    ax.scatter(d.down[~pw], d.up[~pw], s=30, color=ORANGE, edgecolor="white", lw=0.6, zorder=3)
    ax.set_xlim(0.15, 1.0)
    ax.set_ylim(0.15, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("held-out accuracy, decrease sites (AUROC)")
    pname = "decrease" if prim == "down" else "increase"
    ax.text(0.0, 1.02, f"{title}  ·  {int(pw.sum())} of {len(d)} favour {pname}s",
            transform=ax.transAxes, fontsize=8, color=MUTED, ha="left", va="bottom")
    # label a few well-known footprints
    lab = {"KI_MEK": "MEK", "KI_PI3K_mTOR": "PI3K/mTOR", "EGF": "EGF", "TNF": "TNF",
           "DNA_damage": "DNA damage", "KI_AURK_PLK": "Aurora/PLK", "KI_RAF": "RAF", "JAKSTAT_cytokine": "JAK-STAT",
           "starvation": "starvation"}
    for _, r in d.iterrows():
        if r.footprint in lab:
            off = {"KI_mTORcat": (0, -9, "center"), "KI_RAF": (0, 9, "center"), "EGF": (5, -2, "left"), "TNF": (5, -2, "left"),
                   "JAKSTAT_cytokine": (-4, -3, "right"), "starvation": (-4, 3, "right")}.get(r.footprint, (4, 2, "left"))
            ax.annotate(lab[r.footprint], (r.down, r.up), xytext=off[:2], textcoords="offset points", fontsize=6.8,
                        color=MUTED, ha=off[2], va="center")
bI.set_ylabel("held-out accuracy,\nincrease sites (AUROC)")
plt.setp(bS.get_yticklabels(), visible=False)
bI.text(0.97, 0.16, "below the line:\ndecrease sites\ntransfer better", fontsize=6.8, color=MUTED, ha="right",
        va="bottom", transform=bI.transAxes)
bI.text(0.03, 0.97, "above the line:\nincrease sites\ntransfer better", fontsize=6.8, color=MUTED, ha="left",
        va="top", transform=bI.transAxes)

ptitle(fig, aI, "a", "Mirror test: changes in the primary direction replicate across studies", dx=0.09, dy=0.06)
ptitle(fig, bI, "b", "Footprint transfer to held-out studies", dx=0.09, dy=0.06)
h = [Line2D([], [], color=GREEN, marker="o", lw=0, ms=5.5, label="primary direction better"),
     Line2D([], [], color=ORANGE, marker="o", lw=0, ms=5.5, label="reverse direction better")]
fig.legend(handles=h, loc="upper right", bbox_to_anchor=(0.99, 0.505), ncol=2, fontsize=7.3, handletextpad=0.2)
save(fig, "fig5_portability")
