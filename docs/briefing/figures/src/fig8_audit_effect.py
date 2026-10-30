"""Fig 8 - how independent audits and validations changed claim strength."""
import textwrap
from matplotlib.lines import Line2D
from style import plt, POS, CAV, NEG, VAL, INK, GRID, WIDTH, save

# claim label, self score, after score, outcome class, status text
ROWS = [
    ("K1  CK2 → KPNA4 S60", 3, 3.5, "pos",
     "Confidence raised to 'high' by public CK2 rescue data (no new score assigned)"),
    ("K2  Footprint coherence\n      in tumours", 4, 3, "pos",
     '"Qualified" by audit; its ischaemia caveat was later withdrawn by measured ischaemia times (V-B)'),
    ("K3  Portability asymmetry", 4, 3, "cav",
     'Reworded: "not cell-intrinsic"; the inhibitor arm is MEK-dependent'),
    ("K4  Receptor → ERK", 3, 3, "cav",
     'Reworded from "non-transitivity" to "abundance-driven receptor signal is buffered"'),
    ("F   PhosphoSitePlus\n      corrections", 3, 2, "cav", "Downgraded"),
    ("H   Turnover atlas", 3.5, 3, "cav", "Downgraded"),
    ("Immune–phosphatase\nmimicry", 3, 1, "neg", "Failed audit; dropped"),
    ("G   ERK ceiling", 2, 1, "neg",
     "Robust in the discovery cohort; failed replication; dropped"),
]
COL = {"pos": POS, "cav": CAV, "neg": NEG}

n = len(ROWS)
H = 0.55 * n + 1.75
fig = plt.figure(figsize=(WIDTH, H))
bottom_in, top_in = 1.35, 0.50
ax = fig.add_axes([0.27, bottom_in / H, 0.33, 1 - (bottom_in + top_in) / H])

ys = list(range(n))[::-1]
for y, (lab, s0, s1, oc, status) in zip(ys, ROWS):
    c = COL[oc]
    ax.axhline(y, color=GRID, lw=0.8, zorder=0)
    if s0 != s1:
        ax.annotate("", xy=(s1 + (-0.13 if s1 > s0 else 0.13), y), xytext=(s0, y),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8, mutation_scale=10,
                                    shrinkA=5, shrinkB=0), zorder=2)
    ax.plot(s0, y, "o", ms=8, mfc="white", mec="#555555", mew=1.4, zorder=3)
    ax.plot(s1, y, "o", ms=8, mfc=c, mec=c if s0 != s1 else "#555555",
            mew=1.4 if s0 == s1 else 0, zorder=4)
    if s0 == s1:  # same score: show the open ring around the filled dot
        ax.plot(s1, y, "o", ms=11, mfc="none", mec="#555555", mew=1.2, zorder=5)
    ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=8.5, linespacing=1.1,
            fontweight="bold" if lab.startswith("K") else "normal")
    ax.text(1.04, y, textwrap.fill(status, 40), transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=8, color=c if oc != "neg" else "#555555",
            linespacing=1.15)

# K2: caveat raised then withdrawn
ax.text(3.5, ys[1] - 0.33, "caveat raised, then withdrawn", fontsize=7, ha="center",
        va="center", color=POS, fontstyle="italic")

ax.set_xlim(0.6, 5.2)
ax.set_ylim(-0.6, n - 0.4)
ax.set_xticks([1, 2, 3, 4, 5])
ax.set_xlabel('"Compelling" score (1–5)')
ax.set_yticks([])
ax.spines["left"].set_visible(False)
for x in [1, 2, 3, 4, 5]:
    ax.axvline(x, color=GRID, lw=0.6, zorder=0)

fig.text(0.27 + 0.33 + 0.33 * 0.04, 1 - (top_in - 0.06) / H, "Final status", fontsize=9,
         fontweight="bold", va="bottom")
fig.text(0.015, 1 - 0.08 / H, "How independent audits changed claims", fontsize=10, va="top")

handles = [
    Line2D([], [], marker="o", ls="", ms=8, mfc="white", mec="#555555", mew=1.4,
           label="first proposed (self-assessed)"),
    Line2D([], [], marker="o", ls="", ms=8, mfc=POS, mec=POS, label="after audit: upheld or strengthened"),
    Line2D([], [], marker="o", ls="", ms=8, mfc=CAV, mec=CAV, label="narrowed, reworded or downgraded"),
    Line2D([], [], marker="o", ls="", ms=8, mfc=NEG, mec=NEG, label="failed / dropped"),
]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.005, 0.36 / H), ncol=2,
           frameon=False, fontsize=8, handletextpad=0.3, columnspacing=1.1)

fig.text(0.015, 0.012,
         "Non-numeric verdicts plotted as: \"qualified\" = 3; \"fails\" / failed replication = 1.\nK1 self-score was not recorded at V-I03 and is shown as 3.",
         fontsize=7.5, color="#555555", va="bottom", linespacing=1.2)

save(fig, "fig8_audit_effect")
