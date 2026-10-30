"""Fig 1 - the project at a glance: inputs, agent workflow, outcome counts."""
from matplotlib.patches import FancyBboxPatch
from style import plt, POS, CAV, NEG, VAL, INF, HIL, INK, WIDTH, panel_letter, blank, save


def tint(hexc, a=0.15):
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hexc)
    return (1 - a * (1 - r), 1 - a * (1 - g), 1 - a * (1 - b))


def box(ax, x, y, w, h, text, ec, fc=None, fs=8.5, bold=False, ls="-", lw=1.1, ha="center"):
    fc = fc if fc is not None else tint(ec)
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                       fc=fc, ec=ec, lw=lw, ls=ls, mutation_aspect=1)
    ax.add_patch(p)
    tx = x + w / 2 if ha == "center" else x + 0.012
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=INK, linespacing=1.25)
    return p


def arrow(ax, xy0, xy1, color="#555555", rad=0.0, lw=1.1, ls="-"):
    ax.annotate("", xy=xy1, xytext=xy0,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls,
                                shrinkA=0, shrinkB=0, mutation_scale=9,
                                connectionstyle=f"arc3,rad={rad}"))


fig = plt.figure(figsize=(WIDTH, 8.4))
ax_a = fig.add_axes([0.03, 0.665, 0.95, 0.30])
ax_b = fig.add_axes([0.03, 0.255, 0.95, 0.37])
ax_c = fig.add_axes([0.03, 0.02, 0.95, 0.19])
for ax in (ax_a, ax_b, ax_c):
    blank(ax)

# ---------------------------------------------------------------- (a) inputs
ax = ax_a
ax.set_title("Inputs: data and tools available to the agents", loc="left", pad=4)
LW = 0.655                     # width of the data column
ax.text(0.0, 0.95, "Data", fontsize=9, fontweight="bold", color=INK, va="center")
data = [
    ("Perturbation\ncompendium", "1,937 phosphoproteomics perturbation conditions\nfrom 178 studies, incl. 850 kinase-inhibitor conditions", "-"),
    ("CPTAC tumour\ncohorts", "10 tumour cohorts: 1,030 tumours,\n523 normal adjacent tissues", "-"),
    ("Prior knowledge", "PhosphoSitePlus; Kinobeads drug-target\naffinities (Klaeger 2017)", "-"),
    ("Fetched during\nthe project", "DepMap/CCLE · GDSC · CPTAC LUAD confirmatory ·\nChen 2020 LUAD · Menyhart 2023 CK2 rescue · GDC\nischaemia times · Kuster-lab KRAS dose series", (0, (3, 2))),
]
rh, rg = 0.17, 0.022
y = 0.86
for i, (head, body, ls) in enumerate(data):
    hh = rh + (0.05 if i == 3 else 0)
    y -= hh
    box(ax, 0.0, y, LW, hh, "", "#555555", fc="#f4f4f4", ls=ls)
    ax.text(0.012, y + hh / 2, head, fontsize=8.5, fontweight="bold", va="center", linespacing=1.15)
    ax.text(0.19, y + hh / 2, body, fontsize=8, va="center", linespacing=1.25)
    y -= rg

# tools
tx = LW + 0.035
tw = 1 - tx
ax.text(tx, 0.95, "Tools from two earlier\nAI-run projects", fontsize=9,
        fontweight="bold", color=INK, va="center", linespacing=1.1)
tools = [("Kinase–substrate predictor", "covers 503 kinases"),
         ("Kinase-activity estimator", "version core-v3")]
th = 0.33
for j, (head, body) in enumerate(tools):
    yy = 0.80 - th - j * (th + 0.06)
    box(ax, tx, yy, tw, th, "", INF)
    ax.text(tx + 0.015, yy + th * 0.64, head, fontsize=8.5, fontweight="bold", va="center")
    ax.text(tx + 0.015, yy + th * 0.34, body, fontsize=8, va="center")

# ---------------------------------------------------------------- (b) workflow
ax = ax_b
ax.set_title("Agent team and workflow", loc="left", pad=4)
nh = 0.17
r1y, r2y = 0.58, 0.22
# legend for node colours (top of panel)
leg = [(HIL, "orchestrator"), (INF, "infrastructure / planning"), (POS, "discovery"),
       (VAL, "validation / audit"), (INK, "reporting")]
lx = 0.0
for c, t in leg:
    ax.add_patch(FancyBboxPatch((lx, 0.915), 0.018, 0.045, boxstyle="round,pad=0,rounding_size=0.004",
                                fc=tint(c, 0.25) if c != INK else "#f4f4f4", ec=c, lw=1))
    ax.text(lx + 0.025, 0.937, t, fontsize=7.5, va="center")
    lx += 0.025 + 0.0095 * len(t) + 0.035

# row 1 (left to right)
row1 = [("PI agent", HIL, True), ("Literature\nreview", INF, False),
        ("Idea pool\nand ranking", INF, False),
        ("Discovery branches\n(parallel, ≤4 agents\n× 4 threads)", POS, True),
        ("Quick cross-\nvalidation", VAL, False)]
xs1 = [0.0, 0.17, 0.34, 0.51, 0.77]
ws1 = [0.13, 0.13, 0.13, 0.21, 0.155]
for (t, c, b), x, wd in zip(row1, xs1, ws1):
    box(ax, x, r1y, wd, nh if "\n×" not in t else nh + 0.07, t, c, bold=b, fs=8.5)
for i in range(4):
    yc = r1y + nh / 2
    arrow(ax, (xs1[i] + ws1[i] + 0.005, yc), (xs1[i + 1] - 0.005, yc))

# row 2 (right to left)
row2 = [("Independent\naudit", VAL, True), ("Claim\nfreeze", INK, False),
        ("Report", INK, False), ("Fact-check", VAL, False)]
xs2 = [0.77, 0.56, 0.36, 0.16]
ws2 = [0.155, 0.13, 0.13, 0.13]
for (t, c, b), x, wd in zip(row2, xs2, ws2):
    box(ax, x, r2y, wd, nh, t, c, bold=b, fc="#f4f4f4" if c == INK else None)
arrow(ax, (xs1[4] + ws1[4] / 2, r1y - 0.005), (xs2[0] + ws2[0] / 2, r2y + nh + 0.005))
for i in range(3):
    yc = r2y + nh / 2
    arrow(ax, (xs2[i] - 0.005, yc), (xs2[i + 1] + ws2[i + 1] + 0.005, yc))

# feedback: audit narrows or kills claim -> discovery branches
arrow(ax, (xs2[0] + 0.02, r2y + nh + 0.005), (xs1[3] + ws1[3] * 0.55, r1y - 0.005),
      color=CAV, rad=-0.25, lw=1.3)
ax.text(0.638, 0.475, "audit narrows\nor kills claim", color=CAV, fontsize=8,
        ha="right", va="center", fontweight="bold", linespacing=1.15)

# milestone retrospectives -> LESSONS.md -> PI agent
box(ax, 0.0, 0.0, 0.36, 0.12, "Milestone retrospectives → LESSONS.md", INF, fs=8)
arrow(ax, (xs2[3] + ws2[3] * 0.5, r2y - 0.005), (xs2[3] + ws2[3] * 0.5, 0.12 + 0.005), color="#555555")
arrow(ax, (0.04, 0.12 + 0.005), (0.04, r1y - 0.005), color=INF, lw=1.3)
ax.text(0.05, 0.40, "lessons feed\nnext round", color=INF, fontsize=8,
        ha="left", va="center", fontweight="bold", linespacing=1.15)

# ---------------------------------------------------------------- (c) outcomes
ax = ax_c
ax.set_title("Outcome counts", loc="left", pad=4)
ax.text(0.0, 0.92, "Agent tasks run", fontsize=8.5, fontweight="bold", va="center")
stats = [("3", "infrastructure\ntasks", INF), ("18", "discovery\nbranches", POS),
         ("4", "targeted\nvalidations", VAL), ("6", "independent\naudits", VAL),
         ("3", "reviews and\nretrospectives", NEG)]
for i, (n, lab, c) in enumerate(stats):
    x = i * 0.115
    ax.text(x, 0.64, n, fontsize=18, fontweight="bold", color=c, va="center", ha="left")
    ax.text(x, 0.34, lab, fontsize=7.5, va="center", ha="left", linespacing=1.1)

# stacked bar of findings (negative count = negative/dropped tiles in fig 2)
N_NEG = 9
bx0, bw = 0.6, 0.4
segs = [(1, "Tier-1 discovery", HIL), (3, "Tier-2 findings", POS),
        (5, "caveated", CAV), (N_NEG, "negative/dropped", NEG)]
tot = sum(s[0] for s in segs)
ax.text(bx0, 0.92, "Findings", fontsize=8.5, fontweight="bold", va="center")
x = bx0
for n, lab, c in segs:
    wseg = bw * n / tot
    ax.add_patch(plt.Rectangle((x, 0.52), wseg, 0.25, fc=c, ec="white", lw=1.2))
    ax.text(x + wseg / 2, 0.645, str(n), color="white", fontsize=9, fontweight="bold",
            ha="center", va="center")
    x += wseg
for k, (n, lab, c) in enumerate(segs):
    lx = bx0 + (k % 2) * 0.2
    ly = 0.33 - (k // 2) * 0.2
    ax.add_patch(plt.Rectangle((lx, ly - 0.05), 0.014, 0.1, fc=c, ec="none"))
    ax.text(lx + 0.022, ly, f"{n} {lab}", fontsize=8, va="center")

for a, l in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
    panel_letter(fig, a, l, dx=-0.025, dy=0.004)

save(fig, "fig1_overview")
