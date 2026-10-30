"""Shared house style for briefing data figures (fig3-fig7)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREEN = "#1b9e77"   # positive / primary
ORANGE = "#d95f02"  # caveated / second series
GREY = "#7f7f7f"    # negative / null / random
PURPLE = "#7570b3"  # audit / validation
PINK = "#e7298a"    # highlight
GRID = "#e6e6e6"
INK = "#222222"
MUTED = "#555555"

OUT = "/workspace/docs/briefing/figures"
R = "/workspace/results"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.8,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "legend.frameon": False,
    "svg.fonttype": "none",
    "axes.titlelocation": "left",
    "axes.titlepad": 6,
})


def letter(fig, ax, s, dx=0.055, dy=0.012):
    """Bold panel letter at the top-left of the panel (left of the y-axis labels).
    Uses figure coordinates, so call after the axes geometry is fixed (no constrained layout)."""
    p = ax.get_position()
    fig.text(p.x0 - dx, p.y1 + dy, s, fontsize=11, fontweight="bold",
             va="bottom", ha="left", color=INK)


def ygrid(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def xgrid(ax):
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=200)
    fig.savefig(f"{OUT}/{name}.svg")
    print("saved", name)


def ptitle(fig, ax, s, title, dx=0.06, dy=0.015, gap=0.028):
    """Bold panel letter plus a 10-pt panel title, both left-aligned at the panel's left edge."""
    p = ax.get_position()
    fig.text(p.x0 - dx, p.y1 + dy, s, fontsize=11, fontweight="bold", va="baseline", ha="left", color=INK)
    if title:
        fig.text(p.x0 - dx + gap, p.y1 + dy, title, fontsize=10, va="baseline", ha="left", color=INK)
