"""Figures for b_I17 -> results/b_I17/fig*.png"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
RES = "/workspace/results/b_I17"
BLUE, ORANGE, AQUA, GRAY, INK = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.dpi": 150})
A = pd.read_csv(f"{RES}/shift_sites.tsv", sep="\t")
X = pd.read_csv(f"{RES}/extra_shift_by_experiment.tsv", sep="\t")
F = pd.read_parquet(f"{RES}/curves_all.parquet")
allo = {"MK2206", "PD0325901", "refametinib", "selumetinib"}

# Fig 1: shift distributions per experiment
dn = A[A.direction == "down"]
order = X.sort_values("median_shift_down").exp_id.tolist()
fig, ax = plt.subplots(figsize=(7.5, 5.2))
for i, ex in enumerate(order):
    d = dn[dn.exp_id == ex]; x = X[X.exp_id == ex].iloc[0]
    col = ORANGE if x.agent in allo else BLUE
    ax.plot([x.window_lo, x.window_hi], [i, i], color="#e4e3df", lw=7, solid_capstyle="butt", zorder=0)
    q = d["shift"].quantile([.1, .25, .5, .75, .9]).values
    ax.plot([q[0], q[4]], [i, i], color=col, lw=1)
    ax.plot([q[1], q[3]], [i, i], color=col, lw=4, solid_capstyle="round")
    ax.plot(q[2], i, "o", color="white", ms=4, mec=col, mew=1.5)
    dd = d[d.ks_direct | d.self]
    ax.plot(dd["shift"], np.full(len(dd), i + 0.28), "|", color=INK, ms=4, alpha=.7)
ax.axvline(0, color=INK, lw=.8, ls="--")
ax.set_yticks(range(len(order)))
ax.set_yticklabels([f"{X.set_index('exp_id').loc[e,'agent']} · {X.set_index('exp_id').loc[e,'cell_line']} · {int(X.set_index('exp_id').loc[e,'time_min'])}′ (n={X.set_index('exp_id').loc[e,'n_down']})" for e in order], fontsize=7)
ax.set_xlabel("ΔpEC50 = pEC50(site) − pKd_app(target)   (>0: site more sensitive than binding)")
ax.set_title("ATP-site inhibitors (blue): sites respond at 3–100x the Kinobeads Kd; allosteric (orange) MEKi at/below Kd.\nGrey band = observable window; box = IQR; ticks = sites on target / curated target substrates", fontsize=7.5, loc="right", color=INK)
fig.tight_layout(); fig.savefig(f"{RES}/fig1_shift_by_experiment.png"); plt.close(fig)

# Fig 2: drug-level pKd vs median pEC50
dd = X.groupby("agent").agg(pKd=("pKd_target", "first"), med=("median_pEC50_down", "median"), p90=("p90_pEC50_down", "median")).reset_index()
fig, ax = plt.subplots(figsize=(4.2, 4))
for _, r in dd.iterrows():
    c = ORANGE if r.agent in allo else BLUE
    ax.plot(r.pKd, r.med, "o", color=c, ms=7, mec="white", mew=1.5)
    ax.annotate(r.agent, (r.pKd, r.med), xytext=(4, 3), textcoords="offset points", fontsize=7, color=INK)
lim = [4, 9.2]; ax.plot(lim, lim, color=GRAY, lw=.8, ls="--"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("pKd_app of designated target (Klaeger 2017)"); ax.set_ylabel("median pEC50 of down-regulated sites")
ax.set_title("Target Kd does not predict where the\nphosphoproteome responds (ρ=0.26, shuffle p=0.21)", fontsize=8, loc="left", color=INK)
fig.tight_layout(); fig.savefig(f"{RES}/fig2_pKd_vs_pEC50.png"); plt.close(fig)

# Fig 3: Hill calibration + replicate
cal = pd.read_parquet(f"{RES}/hill_calibration.parquet")
rep = pd.read_csv(f"{RES}/extra_dasatinib_A549_replicate_hill.tsv", sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(9, 3))
bins = np.log2([0.25, 0.35, 0.45, 0.55, 0.65, 0.78, 0.92, 1.1, 1.3, 1.55, 1.85, 2.25, 2.75, 3.25, 3.75, 4.5, 5.5])
axs[0].hist(np.log2(A.hill), bins=bins, density=True, color=BLUE, alpha=.85, label="observed (sig. curves)")
axs[0].hist(np.log2(cal.hill_sim_h1), bins=bins, density=True, histtype="step", color=ORANGE, lw=2, label="simulated, true Hill = 1")
axs[0].set_xlabel("log2 Hill slope"); axs[0].set_ylabel("density"); axs[0].legend(fontsize=7, frameon=False)
axs[0].set_title("Excess steep curves vs Hill=1 null", fontsize=8, loc="left")
axs[1].scatter(rep.pEC50_a, rep.pEC50_b, s=10, color=BLUE, alpha=.7, edgecolor="none")
axs[1].plot([5, 9], [5, 9], color=GRAY, lw=.8, ls="--"); axs[1].set_xlabel("pEC50 replicate 1"); axs[1].set_ylabel("pEC50 replicate 2")
axs[1].set_title("Dasatinib A549 replicates: ρ=0.75", fontsize=8, loc="left")
j = np.random.default_rng(0).normal(0, .06, (2, len(rep)))
axs[2].scatter(np.log2(rep.hill_a) + j[0], np.log2(rep.hill_b) + j[1], s=10, color=BLUE, alpha=.7, edgecolor="none")
axs[2].set_xlabel("log2 Hill replicate 1"); axs[2].set_ylabel("log2 Hill replicate 2")
axs[2].set_title("Hill: ρ=0.44 (6 vs 1.6 exp. steep in both)", fontsize=8, loc="left")
fig.tight_layout(); fig.savefig(f"{RES}/fig3_hill.png"); plt.close(fig)

# Fig 4: K562 BCR-ABL amplification
S = pd.read_csv(f"{RES}/amplification_sites.tsv", sep="\t")
fig, axs = plt.subplots(1, 2, figsize=(8, 3.6), sharex=True)
for ax, ex, t in zip(axs, ["DS0040_Dasatinib0nM_0", "DS0041_Imatinib0nM_0"], ["dasatinib, K562", "imatinib, K562"]):
    d = S[S.exp_id == ex].copy(); d["grp"] = np.where(d.self2, "BCR-ABL / target", "downstream readout")
    d = d.sort_values(["grp", "pEC50"])
    for i, (_, r) in enumerate(d.iterrows()):
        c = GRAY if r.self2 else AQUA
        ax.plot(r.pEC50, i, "o", color=c, ms=5)
        ax.text(9.6, i, f"{r.gene} {r.site.split('_')[1]}", fontsize=6, va="center", color=INK)
    for g, c in [(True, GRAY), (False, AQUA)]:
        ax.axvline(d.loc[d.self2 == g, "pEC50"].median(), color=c, lw=1.2, ls="--")
    ax.axvline(d.pKd_target.iloc[0], color=ORANGE, lw=1, ls=":")
    ax.set_yticks([]); ax.set_xlim(5, 10.6); ax.set_xlabel("pEC50"); ax.set_title(t + "  (grey=target sites, green=readouts, orange=Kd)", fontsize=8, loc="left")
fig.tight_layout(); fig.savefig(f"{RES}/fig4_bcrabl_amplification.png"); plt.close(fig)
print("ok")
