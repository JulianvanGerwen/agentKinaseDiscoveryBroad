"""fig6: Finding K4 - offsets not gains; buffered receptor abundance. Re-plots saved tables only."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from style import *

t2 = pd.read_csv(f"{R}/b_I24/tables/T2_gain_meta.tsv", sep="\t")
t2 = t2[t2.variant == "tc_padj"].set_index("edge")
t4 = pd.read_csv(f"{R}/b_I24/tables/T4_modifiers_meta.tsv", sep="\t")
t4 = t4[(t4.model == "main") & (t4.primary == True)]
t6 = pd.read_csv(f"{R}/b_I24/tables/T6_spikein_power.tsv", sep="\t")
ej = pd.read_csv(f"{R}/val_E/E_joint_by_cohort.tsv", sep="\t")
ej = ej[ej.spec == "base (purity, prolif)"]
em = pd.read_csv(f"{R}/val_E/E_joint_meta.tsv", sep="\t").set_index("spec")

fig = plt.figure(figsize=(7.2, 6.0))
aA1 = fig.add_axes([0.06, 0.69, 0.13, 0.20])
aA2 = fig.add_axes([0.235, 0.69, 0.13, 0.20])
aB = fig.add_axes([0.625, 0.645, 0.27, 0.285])
aC1 = fig.add_axes([0.255, 0.115, 0.2, 0.33])
aC2 = fig.add_axes([0.475, 0.115, 0.2, 0.33], sharey=aC1)
aD = fig.add_axes([0.80, 0.115, 0.18, 0.33])

# ---------- (a) schematic ----------
x = np.linspace(0, 1, 20)
for ax, kind in [(aA1, "offset"), (aA2, "gain")]:
    ax.plot(x, 0.2 + 0.5 * x, color=GREY, lw=1.6, ls=(0, (3, 2)))
    if kind == "offset":
        ax.plot(x, 0.45 + 0.5 * x, color=PINK, lw=2.2)
        ax.annotate("", xy=(0.5, 0.70), xytext=(0.5, 0.45), arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
        ax.set_title("offset change", fontsize=9, loc="center", pad=3)
        ax.text(0.5, -0.2, "same slope,\nline shifted up", transform=ax.transAxes, ha="center", va="top",
                fontsize=7.3, color=INK)
    else:
        ax.plot(x, 0.2 + 0.95 * x, color=PINK, lw=2.2)
        ax.set_title("gain change", fontsize=9, loc="center", pad=3)
        ax.text(0.5, -0.2, "steeper slope:\nmore output per input", transform=ax.transAxes, ha="center",
                va="top", fontsize=7.3, color=INK)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("upstream input", fontsize=7.5, labelpad=2)
    ax.set_ylabel("downstream output", fontsize=7.5, labelpad=2)
fig.text(0.03, 0.575, "dashed: reference tumours;  pink: tumours with a lesion", fontsize=7, color=MUTED)

# ---------- (b) gains per link ----------
ax = aB
xgrid(ax)
EDGES = [("RTK>ADAPT", "receptor pY → adaptor pY"), ("ADAPT>ERK", "adaptor pY → ERK"),
         ("ERK>RSK", "ERK → RSK"), ("ERK>ERKSUB", "ERK → ERK substrates"),
         ("ADAPT>AKTSUB", "adaptor pY → AKT substrates"), ("AKTSUB>MTORC1", "AKT substrates → mTORC1"),
         ("MTORC1>S6", "mTORC1 → S6"), ("RTK>ERK", "receptor pY → ERK (direct)")]
y = np.arange(len(EDGES))[::-1].astype(float)
y[-1] -= 0.4
for yy, (e, lab) in zip(y, EDGES):
    r = t2.loc[e]
    col = ORANGE if e == "RTK>ERK" else GREEN
    ax.plot([r.std_lo, r.std_hi], [yy, yy], color=col, lw=1.6)
    ax.scatter(r.std_slope_meta, yy, s=30, color=col, edgecolor="white", lw=0.6, zorder=3)
    ax.text(0.665, yy, f"{r.std_slope_meta:.2f}  ({int(r.k)})", va="center", fontsize=7, color=INK)
prod = t2.loc["RTK>ADAPT", "std_slope_meta"] * t2.loc["ADAPT>ERK", "std_slope_meta"]
ax.scatter(prod, y[-1], marker="|", s=110, color=INK, lw=1.4, zorder=4)
ax.annotate(f"expected if passed on\n(product of 2 links) {prod:.2f}", (prod, y[-1]), xytext=(prod + 0.03, y[-1] - 1.05),
            fontsize=6.6, color=INK, va="center", arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
ax.axvline(0, color="#999999", lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([l for _, l in EDGES], fontsize=7.5)
ax.set_xlim(-0.08, 0.66)
ax.set_ylim(y[-1] - 1.7, y[0] + 0.6)
ax.set_xlabel("gain per link (standardised slope)")
ax.text(0.665, y[0] + 0.75, "est. (cohorts)", fontsize=6.6, color=MUTED, va="bottom")

# ---------- (c) lesions: gain change vs offset ----------
ROWS = [("RTK>ERK", "KRAS", "mut", "KRAS mutation", "receptor → ERK"),
        ("ADAPT>AKTSUB", "PTEN", "mut", "PTEN mutation", "adaptor → AKT subs."),
        ("ADAPT>AKTSUB", "PTEN", "cnv", "PTEN copy number", "adaptor → AKT subs."),
        ("AKTSUB>MTORC1", "TSC1", "protein", "TSC1 protein", "AKT subs. → mTORC1"),
        ("AKTSUB>MTORC1", "TSC2", "protein", "TSC2 protein", "AKT subs. → mTORC1"),
        ("ADAPT>ERK", "DUSP6", "protein", "DUSP6 protein", "adaptor → ERK"),
        ("ADAPT>ERK", "NF1", "mut", "NF1 mutation", "adaptor → ERK"),
        ("ADAPT>ERK", "DUSP4", "cnv", "DUSP4 copy number", "adaptor → ERK")]
yc = np.arange(len(ROWS))[::-1]
pw = t6[(t6.power >= 0.8)].groupby("gene").delta.min()
band = pw.min()
for ax in (aC1, aC2):
    xgrid(ax)
    ax.axvline(0, color="#999999", lw=0.8)
aC1.add_patch(Rectangle((-band, -0.8), 2 * band, len(ROWS) + 0.3, facecolor=PURPLE, alpha=0.13, edgecolor="none",
                        zorder=0))
aC1.text(0, len(ROWS) - 0.35, f"±{band:.2f}: detectable\nwith ≥84% power", ha="center", va="bottom", fontsize=6.6,
         color=PURPLE)
for yy, (e, g, ty, lab, edge) in zip(yc, ROWS):
    r = t4[(t4.edge == e) & (t4.gene == g) & (t4.type == ty)].iloc[0]
    aC1.plot([r.int_lo, r.int_hi], [yy, yy], color=GREY, lw=1.5)
    aC1.scatter(r.int_meta, yy, s=26, color=GREY, edgecolor="white", lw=0.6, zorder=3)
    sig = r.lev_q_primary < 0.05
    col = GREEN if sig else GREY
    aC2.plot([r.lev_lo, r.lev_hi], [yy, yy], color=col, lw=1.5)
    aC2.scatter(r.lev_meta, yy, s=26, color=col, edgecolor="white", lw=0.6, zorder=3)
    q = r.lev_q_primary
    qs = f"q {q:.0e}".replace("e-0", "e-") if q < 0.001 else (f"q {q:.3f}" if q < 0.01 else f"q {q:.2f}")
    aC2.text(1.16, yy, qs, fontsize=6.5, va="center",
             color=INK if sig else MUTED)
aC1.set_yticks(yc)
aC1.set_yticklabels([f"{lab}\n({edge})" for _, _, _, lab, edge in ROWS], fontsize=7, linespacing=0.95)
plt.setp(aC2.get_yticklabels(), visible=False)
aC1.set_ylim(-0.7, len(ROWS) + 0.6)
aC1.set_xlim(-0.75, 0.75)
aC2.set_xlim(-0.75, 1.13)
aC1.set_xticks([-0.5, 0, 0.5])
aC2.set_xticks([-0.5, 0, 0.5, 1.0])
aC1.set_xlabel("change in gain\n(slope difference)")
aC2.set_xlabel("change in offset\n(SD of downstream output)")
aC1.text(0.5, 1.01, "gain", transform=aC1.transAxes, ha="center", va="bottom", fontsize=8.5, fontweight="bold",
         color=INK)
aC2.text(0.5, 1.01, "offset", transform=aC2.transAxes, ha="center", va="bottom", fontsize=8.5,
         fontweight="bold", color=INK)
ngain = (t4.int_q_primary < 0.05).sum()
fig.text(0.255, 0.477, f"{ngain} of {len(t4)} pre-declared gain tests significant (smallest q {t4.int_q_primary.min():.2f})",
         fontsize=7.2, color=MUTED)

# ---------- (d) receptor abundance ----------
ax = aD
xgrid(ax)
order = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD"]
ej = ej.set_index("cohort").loc[order]
yd = np.arange(len(order))[::-1] + 1.6
for yy, (c, r) in zip(yd, ej.iterrows()):
    lo, hi = r.b_prot - 1.96 * r.se_prot, r.b_prot + 1.96 * r.se_prot
    ax.plot([lo, hi], [yy, yy], color=GREY, lw=1.3)
    ax.scatter(r.b_prot, yy, s=18 + r.n / 6, color=GREY, edgecolor="white", lw=0.5, zorder=3)
mb = em.loc["base (purity, prolif)"]
ax.plot([mb.prot_lo, mb.prot_hi], [0.6, 0.6], color=GREEN, lw=2.2)
ax.scatter(mb.prot, 0.6, marker="D", s=46, color=GREEN, edgecolor="white", lw=0.6, zorder=3)
ax.text(mb.prot, 0.0, f"{mb.prot:.2f} [{mb.prot_lo:.2f}, {mb.prot_hi:.2f}]".replace("-", "−") + f"\np = {mb.prot_p:.3f}",
        ha="center", va="top", fontsize=6.8, color=INK)
ax.axvline(0, color="#999999", lw=0.8)
ax.set_yticks(list(yd) + [0.6])
ax.set_yticklabels(order + ["meta (7)"], fontsize=7.5)
ax.set_ylim(-1.2, yd[0] + 0.7)
ax.set_xlim(-0.7, 0.35)
ax.set_xticks([-0.6, -0.3, 0, 0.3])
ax.set_xlabel("ERK activation per SD of\nreceptor protein (SD)")

ptitle(fig, aA1, "a", "Two ways a lesion can act", dx=0.05, dy=0.045)
ptitle(fig, aB, "b", "Each link of the cascade transmits", dx=0.215, dy=0.03)
ptitle(fig, aC1, "c", "Lesions shift offsets, not gains", dx=0.245, dy=0.06)
ptitle(fig, aD, "d", "More receptor, less ERK", dx=0.065, dy=0.06)
fig.text(0.735, 0.477, "at fixed receptor pY", fontsize=7.2, color=MUTED)
save(fig, "fig6_cascade")
