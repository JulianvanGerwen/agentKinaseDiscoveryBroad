import os
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS"]: os.environ[v]="4"
import sys; sys.path.insert(0, "/workspace/scripts/b_I23")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from importlib import import_module
P = import_module("04_primary")
T = "/workspace/results/b_I23/tables/"; FG = "/workspace/results/b_I23/figures/"
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
D = P.D.reset_index(drop=True)
Y = P.resid(D, "auc"); X = P.resid(D, "LR")
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1, 1.15]})
a = ax[0]; a.axhline(0, color="#ddd", lw=1); a.axvline(0, color="#ddd", lw=1)
a.scatter(X, Y, s=28, color=BLUE, edgecolor="#fcfcfb", lw=1.5, zorder=3)
b = np.polyfit(X, Y, 1); xx = np.linspace(X.min(), X.max(), 10); a.plot(xx, np.polyval(b, xx), color=INK, lw=2)
top = np.argsort(-np.abs(Y))[:7]
for i in top: a.annotate(f"{D.drug[i]}·{D.cell[i]}", (X[i], Y[i]), fontsize=7, color=MUTED, xytext=(4, 2), textcoords="offset points")
a.set_xlabel("reverse/primary log-ratio LR (drug & cell FE removed)"); a.set_ylabel("GDSC AUC (drug & cell FE removed)\n← more sensitive | more resistant →")
a.set_title("P1 (Hijazi, 32 drugs × 3 lines): std β −0.28 [−0.47, −0.05], p=0.025", fontsize=9, color=INK, loc="left")
R4 = pd.read_csv(T + "T4_primary_results.tsv", sep="\t"); R5 = pd.read_csv(T + "T5_robustness.tsv", sep="\t"); R8 = pd.read_csv(T + "T8_baseline_compare.tsv", sep="\t")
K = pd.read_csv(T + "T6_kinase_groups.tsv", sep="\t"); K = K[K.feature == "z"]
rows = [("P1 primary (LR → AUC)", R4.iloc[0]), ("LR → ln IC50", R4.iloc[1]), ("log f_up | f_down", R4.iloc[2]), ("log f_down | f_up", R4.iloc[3]),
        ("LR | target engagement", R4.iloc[7]), ("LR | DepMap target expr", R8.iloc[1])]
rows += [(f"{lab}", R5[R5.test == t].iloc[0]) for lab, t in [("GDSC2 only", "GDSC2_only"), ("GDSC1 only", "GDSC1_only"), ("drop HL60", "drop_HL60"), ("drop MCF7", "drop_MCF7"), ("drop NTERA2", "drop_NTERA2")]]
rows += [(f"z {r.group} (S3)", r) for _, r in K.iterrows()]
a = ax[1]; a.axvline(0, color="#bbb", lw=1)
for j, (lab, r) in enumerate(rows):
    c = ORANGE if j == 0 else (BLUE if "S3" not in lab else MUTED)
    a.plot([r.ci_lo, r.ci_hi], [j, j], color=c, lw=2); a.scatter([r.std_beta], [j], color=c, s=30, zorder=3)
a.set_yticks(range(len(rows))); a.set_yticklabels([l for l, _ in rows], fontsize=7.5); a.invert_yaxis()
a.set_xlabel("standardized β on AUC (drug-cluster bootstrap 95% CI)\n< 0: larger feature ↔ more sensitive")
a.set_title("Robustness and kinase-group tests", fontsize=9, color=INK, loc="left")
fig.tight_layout(); fig.savefig(FG + "fig1_reverse_layer_vs_sensitivity.png", dpi=160)
print("ok")
