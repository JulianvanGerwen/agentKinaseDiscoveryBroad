"""Fig 2 - research map: every branch / validation / audit and its outcome."""
from matplotlib.patches import FancyBboxPatch, Rectangle
from style import plt, POS, CAV, NEG, VAL, INF, HIL, INK, WIDTH, save
import matplotlib.colors as mc


def tint(hexc, a=0.16):
    r, g, b = mc.to_rgb(hexc)
    return (1 - a * (1 - r), 1 - a * (1 - g), 1 - a * (1 - b))


OUTCOME = {  # colour, glyph, legend label
    "pos": (POS, "✓", "positive → key finding (K1–K4)"),
    "cav": (CAV, "~", "caveated"),
    "neg": (NEG, "✗", "negative / dropped"),
    "val": (VAL, "◆", "validation / audit"),
}

# theme -> rows -> tiles (id, name, outcome, note, span)
THEMES = [
    ("Inhibitor and\nstimulus\nresponses", [
        [("I19", "Feedback atlas", "pos", "→ K3", 1),
         ("A-I19", "Audit of feedback atlas", "val", "", 1),
         ("I05", "Timing of responses", "neg", "supports K3", 1)],
        [("I22", "Reverse layer\ncell-intrinsic?", "neg", "", 1),
         ("I23", "Reverse layer vs\ndrug resistance", "neg", "", 1),
         ("I25", "Reverse-layer\nmechanism", "neg", "", 1)],
    ]),
    ("Dose\npharmacology", [
        [("I17", "Dose vs drug affinity", "cav", "ATP consistency", 1)],
    ]),
    ("Kinase–substrate\nrelationships", [
        [("I03", "Integrated new edges;\nmTORC1 candidates", "pos", "→ K1", 1),
         ("V-I03", "Discriminating\nvalidation", "val", "", 1),
         ("A-D", "Audit of new edges", "val", "", 1)],
        [("I26", "PhosphoSitePlus audit", "cav", "", 1),
         ("A-F", "Audit + public\nCK2 rescue", "val", "", 1),
         ("I28", "Turnover atlas", "cav", "", 1),
         ("A-H", "Audit of\nturnover atlas", "val", "", 1)],
        [("I01", "Dark kinases", "neg", "", 1),
         ("I20", "Copy-number\ninstruments", "neg", "", 1)],
    ]),
    ("Tumours", [
        [("I10", "Perturbation footprints\nin tumours", "pos", "→ K2", 1),
         ("A-I10", "Audit of footprints", "val", "", 1),
         ("V-B", "Measured ischaemia", "val", "", 1),
         ("Q8", "Tumour vs normal\ntissue", "pos", "→ K2", 1)],
        [("I24", "Cascade gain", "pos", "→ K4", 1),
         ("A-I24", "Audit of cascade gain", "val", "", 1),
         ("V-E", "Ischaemia re-check", "val", "", 1),
         ("I08", "Abundance vs\nactivating sites", "cav", "partly known", 1)],
        [("I27", "Allele-specific\nsignalling", "neg", "", 1),
         ("I29", "ERK ceiling in\nEGFR-mutant LUAD", "neg", "failed replication", 1),
         ("I13", "Phosphatase lesions", "neg", "", 1)],
    ]),
    ("Single-dataset\ncandidates", [
        [("V1", "Quick cross-validation", "val",
          "killed MET/ALK → Hippo, PAK4 S148 and most reporter sites", 4)],
    ]),
]

# ------------------------------------------------------------------ geometry (inches)
LABW = 1.42
NCOL = 4
GAP = 0.06
TILEW = (WIDTH - LABW - 0.05 - (NCOL - 1) * GAP) / NCOL
TILEH = 0.50
RGAP = 0.07
THGAP = 0.20
LEGH = 0.40
TLH = 1.35

nrows = sum(len(r) for _, r in THEMES)
mapH = nrows * TILEH + (nrows - len(THEMES)) * RGAP + (len(THEMES) - 1) * THGAP
H = 0.35 + LEGH + mapH + 0.35 + TLH
fig = plt.figure(figsize=(WIDTH, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, WIDTH)
ax.set_ylim(0, H)
ax.set_axis_off()

top = H - 0.08
ax.text(0.05, top, "a", fontsize=11, fontweight="bold", va="top")
ax.text(0.25, top, "Every line of work and its outcome", fontsize=10, va="top")

# legend
ly = top - 0.42
lx = LABW
for key in ("pos", "cav", "neg", "val"):
    c, g, lab = OUTCOME[key]
    ax.add_patch(FancyBboxPatch((lx, ly - 0.09), 0.18, 0.18,
                                boxstyle="round,pad=0,rounding_size=0.03",
                                fc=tint(c), ec=c, lw=1.2))
    ax.text(lx + 0.09, ly, g, color=c, fontsize=9, fontweight="bold", ha="center", va="center")
    ax.text(lx + 0.25, ly, lab, fontsize=8, va="center")
    lx += 0.25 + 0.06 * len(lab) + 0.15

counts = {k: 0 for k in OUTCOME}
y = ly - 0.30
for ti, (theme, rows) in enumerate(THEMES):
    ytop = y
    for r, row in enumerate(rows):
        x = LABW
        for (tid, name, oc, note, span) in row:
            counts[oc] += 1
            c, g, _ = OUTCOME[oc]
            w = TILEW * span + GAP * (span - 1)
            ax.add_patch(FancyBboxPatch((x, y - TILEH), w, TILEH,
                                        boxstyle="round,pad=0,rounding_size=0.05",
                                        fc=tint(c), ec=c, lw=1.2))
            ax.add_patch(Rectangle((x, y - TILEH + 0.03), 0.05, TILEH - 0.06, fc=c, ec="none"))
            ax.text(x + 0.11, y - 0.06, tid, fontsize=8.5, fontweight="bold", va="top", color=INK)
            ax.text(x + w - 0.07, y - TILEH + 0.05, g, fontsize=9.5, fontweight="bold",
                    va="bottom", ha="right", color=c)
            if span > 1:
                ax.text(x + 0.11, y - 0.24, name, fontsize=8, va="top")
                ax.text(x + 1.75, y - 0.24, note, fontsize=7.5, va="top", style="italic", color=INK)
            else:
                nm_lines = name.count("\n") + 1
                ax.text(x + 0.40 if nm_lines == 1 and False else x + 0.11,
                        y - 0.205, name, fontsize=7.5, va="top", linespacing=1.1)
                if note:
                    # note sits right of the ID on the first line
                    ax.text(x + 0.11 + 0.075 * len(tid) + 0.08, y - 0.07, note, fontsize=7,
                            va="top", style="italic", color=c if oc != "neg" else "#555555",
                            fontweight="bold" if note.startswith("→") else "normal")
            x += w + GAP
        y -= TILEH + (RGAP if r < len(rows) - 1 else 0)
    # theme label + bracket
    ax.plot([LABW - 0.08, LABW - 0.08], [y + 0.02, ytop - 0.02], color="#999999", lw=1)
    ax.text(LABW - 0.15, (y + ytop) / 2, theme, fontsize=8.5, fontweight="bold",
            ha="right", va="center", linespacing=1.15)
    y -= THGAP

print("tile counts:", counts)

# ------------------------------------------------------------------ timeline (b)
ty0 = y + THGAP - 0.35
ax.text(0.05, ty0, "b", fontsize=11, fontweight="bold", va="top")
ax.text(0.25, ty0, "Timeline (times approximate, hours since kickoff)", fontsize=10, va="top")
x0, x1 = LABW, WIDTH - 0.35
tmax = 18
lineY = ty0 - 0.62


def tx(t):
    return x0 + (x1 - x0) * t / tmax


ax.plot([x0, x1], [lineY, lineY], color="#bbbbbb", lw=2.5, solid_capstyle="round", zorder=1)
for h in range(0, tmax + 1, 2):
    ax.plot([tx(h), tx(h)], [lineY - 0.04, lineY + 0.04], color="#bbbbbb", lw=0.8, zorder=1)
events = [(0, "T+0", "kickoff and\nliterature", INF, "below"),
          (1, "≈1 h", "infrastructure ready", INF, "above"),
          (6, "≈6 h", "milestone M1 retro", NEG, "above"),
          (14, "≈14 h", "milestone M2 retro", NEG, "above"),
          (18, "≈18 h", "claims frozen\nand report", HIL, "below")]
for t, lab, desc, c, pos in events:
    ax.plot(tx(t), lineY, "o", ms=7, mfc=c, mec="white", mew=1, zorder=3)
    if pos == "above":
        ax.text(tx(t), lineY + 0.12, f"{lab} {desc}", fontsize=7.5, ha="left" if t < 17 else "right",
                va="bottom", color=INK)
    else:
        ax.text(tx(t), lineY - 0.13, f"{lab}\n{desc}", fontsize=7.5,
                ha="left" if t < 1 else "right", va="top", color=INK, linespacing=1.1)
# branches-running band
ax.add_patch(Rectangle((tx(1), lineY - 0.025), tx(18) - tx(1), 0.05, fc=POS, alpha=0.55, ec="none", zorder=2))
ax.text(tx(9.5), lineY - 0.13, "discovery branches, audits and validations run in parallel",
        fontsize=7.5, ha="center", va="top", color=POS)

save(fig, "fig2_research_map")
