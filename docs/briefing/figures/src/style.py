"""Shared briefing figure style (process figures)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

POS = "#1b9e77"      # key finding / positive
CAV = "#d95f02"      # caveated / partial
NEG = "#7f7f7f"      # negative / dropped
VAL = "#7570b3"      # validation / audit
INF = "#66a61e"      # infrastructure
HIL = "#e7298a"      # highlight
GRID = "#e6e6e6"
INK = "#222222"
WIDTH = 7.2

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#444444",
    "axes.grid": False,
    "grid.color": GRID,
    "svg.fonttype": "none",
    "text.color": INK,
    "axes.labelcolor": INK,
})


def panel_letter(fig, ax, letter, dx=-0.0, dy=0.0):
    """Bold panel letter at top-left of the axes bbox (figure coords)."""
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=11, fontweight="bold",
             ha="left", va="bottom")


def blank(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=200)
    print("saved", name)
