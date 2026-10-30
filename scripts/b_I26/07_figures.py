"""I26 figures: fig1 audit vs random control (+ sensitivity), fig2 PSP evidence strata, fig3 activity benchmark."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *

BLUE, ORANGE, GRAY, INK, INK2 = "#2a78d6", "#eb6834", "#9a9a94", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
T = f"{R26}/tables"

# ---- fig1
pg = pd.read_csv(f"{T}/per_group_support.tsv", sep="\t", index_col=0)
pg = pg[pg.n_testable >= 15].sort_values("n_testable", ascending=False)
G = pd.read_csv(f"{T}/sensitivity_grid.tsv", sep="\t")
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1.6]))
x = np.arange(len(pg)); w = 0.38
for a, col, title in [(ax[0], "support_rate", "A  Supported (consistent decrease)"), (ax[1], "unsupported_rate", "B  Powered-unsupported")]:
    a.bar(x - w / 2 - 0.01, pg[col], w, color=BLUE, label="PSP edges")
    a.bar(x + w / 2 + 0.01, pg["rand_" + col], w, color=GRAY, label="random sites (same group, matched)")
    a.set_xticks(x, [f"{g}\n(n={n})" for g, n in zip(pg.index, pg.n_testable)], fontsize=7, rotation=30, ha="right")
    a.set_ylabel("fraction of testable edges"); a.set_title(title, loc="left", fontsize=10)
ax[0].legend(frameon=False, fontsize=8)
G["lab"] = "ctrl " + G.ctrl_thr.astype(str) + " / dec " + G.dec_thr.astype(str)
g3 = G[G.gold_min_LT == 3].reset_index(drop=True)
y = np.arange(len(g3))
for c, mk in [("supported", "o"), ("unsupported", "s")]:
    ax[2].scatter(g3[f"psp_{c}"], y, color=BLUE, marker=mk, s=30, label=f"PSP {c}", zorder=3)
    ax[2].scatter(g3[f"rand_{c}"], y, color=GRAY, marker=mk, s=30, label=f"random {c}", zorder=3)
    for i in y:
        ax[2].plot([g3[f"psp_{c}"][i], g3[f"rand_{c}"][i]], [i, i], color="#d0cfca", lw=1, zorder=1)
ax[2].set_yticks(y, [f"{l}  (n={n})" for l, n in zip(g3.lab, g3.n_testable)], fontsize=7)
ax[2].set_xlabel("fraction of testable edges"); ax[2].set_title("C  Threshold sensitivity", loc="left", fontsize=10)
ax[2].legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.tight_layout(); fig.savefig(f"{R26}/figures/fig1_audit_vs_random.png", dpi=160); plt.close()

# ---- fig2 strata
S = pd.read_csv(f"{T}/relations_strata.tsv", sep="\t")
fig, ax = plt.subplots(1, 3, figsize=(11, 3.2), sharey=True)
for a, var, ttl in zip(ax, ["LT_bin", "motif_bin", "evidence"], ["A  PSP low-throughput refs (site LT_LIT)", "B  Kinase motif percentile (Cantley)", "C  PSP edge evidence type"]):
    d = S[S.variable == var]
    xx = np.arange(len(d))
    a.bar(xx - 0.2, d.support_rate, 0.38, color=BLUE, label="supported")
    a.bar(xx + 0.2, d.unsupported_rate, 0.38, color=ORANGE, label="unsupported")
    a.set_xticks(xx, [f"{l}\n(n={n})" for l, n in zip(d.level, d.n)], fontsize=8); a.set_title(ttl, loc="left", fontsize=9.5)
ax[0].set_ylabel("fraction of testable (group, site) units"); ax[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{R26}/figures/fig2_evidence_strata.png", dpi=160); plt.close()

# ---- fig3 activity
A = pd.read_csv(f"{T}/activity_benchmark_summary.tsv", sep="\t")
A = A[(A.targets == "high+medium")]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.2), sharey=True)
vs = ["prune", "prune+amb", "prune+add", "rand", "rand_any"]
labs = ["drop unsupported", "drop unsup.+ambiguous", "drop unsup. + add corrections", "random drop (non-supported pool)", "random drop (any edge)"]
for a, sub in zip(ax, ["all", "affected"]):
    for k, (sch, col) in enumerate([("LOSO", BLUE), ("half-split x20", ORANGE)]):
        d = A[(A.scheme == sch) & (A.subset == sub)].set_index("variant").reindex(vs)
        a.barh(np.arange(len(vs)) + (k - 0.5) * 0.38, d.mean_delta, 0.36, color=col, label=sch)
    a.axvline(0, color=INK2, lw=0.8); a.set_yticks(np.arange(len(vs)), labs)
    n = A[(A.scheme == "LOSO") & (A.subset == sub)].n_cond_evals.iloc[0]
    a.set_xlabel("mean Δ AUROC vs full PSP library"); a.set_title(f"{'A' if sub == 'all' else 'B'}  {sub} inhibitor conditions (LOSO n={n})", loc="left", fontsize=9.5)
ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{R26}/figures/fig3_activity_benchmark.png", dpi=160); plt.close()
print("ok")
