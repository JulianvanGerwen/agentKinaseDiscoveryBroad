"""fig7: three caveated observations. Re-plots saved tables only."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from style import *

ce = pd.read_csv(f"{R}/audit_D/tables/class_effects_loso.tsv", sep="\t")
h1 = pd.read_csv(f"{R}/audit_H/h1_lab_splits.tsv", sep="\t")
nm = pd.read_csv(f"{R}/audit_H/h1_noise_matched_random.tsv", sep="\t")
atp = pd.read_csv(f"{R}/b_I17/atp_drug_table.tsv", sep="\t")
imp = pd.read_csv(f"{R}/b_I17/atp_implied.tsv", sep="\t").set_index("agent")

fig = plt.figure(figsize=(7.2, 6.0))
aA = fig.add_axes([0.275, 0.525, 0.30, 0.39])
aB = fig.add_axes([0.745, 0.585, 0.24, 0.33])
aC = fig.add_axes([0.24, 0.075, 0.42, 0.235])

# ---------- (a) mTORC1 candidates ----------
ax = aA
xgrid(ax)
SITES = [("LAMTOR1 S98", "cand"), ("LAMTOR1 S97", "cand"), ("RRAGC S15", "cand"),
         ("PRAS40 S183", "direct"), ("4EBP1 S65", "direct"), ("4EBP1 T46", "direct"), ("LARP1 S766", "direct"),
         ("EIF4B S422", "s6k"), ("RPS6 S240", "s6k")]
CLS = [("cat_mTORi", "mTOR kinase inhibitors", GREEN, "o", True),
       ("AKTi", "AKT inhibitors", ORANGE, "s", True),
       ("S6Ki", "S6K inhibitors (inhibition not verified)", GREY, "D", False)]
y = []
yy = 0.0
prev = None
for s, g in SITES:
    if prev is not None and g != prev:
        yy -= 0.7
    y.append(yy); yy -= 1; prev = g
y = np.array(y)
piv = ce.set_index(["site", "cls"])
for yi, (s, g) in zip(y, SITES):
    vals = [piv.loc[(s, c)] if (s, c) in piv.index else None for c, *_ in CLS]
    xs = [v.med_z for v in vals if v is not None]
    ax.plot([min(xs), max(xs)], [yi, yi], color="#cccccc", lw=1, zorder=1)
    for (c, lab, col, mk, filled), v in zip(CLS, vals):
        if v is None:
            continue
        ax.scatter(v.med_z, yi, s=30, marker=mk, color=col if filled else "white", edgecolor=col,
                   lw=1.1 if not filled else 0.5, zorder=3)
    m = piv.loc[(s, "cat_mTORi")]
    ax.text(1.7, yi, f"{int(m.n_units)}", fontsize=6.8, va="center", ha="center", color=MUTED)
ax.text(1.7, y[0] + 0.55, "n", fontsize=6.5, ha="center", va="bottom", color=MUTED, linespacing=0.95)
ax.axvline(0, color="#999999", lw=0.8)
ax.axvline(-2, color=GREY, lw=0.8, ls=(0, (3, 2)))
ax.set_yticks(y)
ax.set_yticklabels([s for s, _ in SITES], fontsize=7.8)
for t, (s, g) in zip(ax.get_yticklabels(), SITES):
    if g == "cand":
        t.set_color(PINK); t.set_fontweight("bold")
ax.set_ylim(y[-1] - 0.7, y[0] + 1.3)
ax.set_xlim(-11.5, 2.2)
ax.set_xlabel("median change in phosphorylation\n(z-score within experiment)")
gl = {"cand": "candidate\nmTORC1 sites", "direct": "known direct\nmTORC1 sites", "s6k": "S6K sites\n(controls)"}
for g in ["cand", "direct", "s6k"]:
    yg = y[[i for i, (_, gg) in enumerate(SITES) if gg == g]]
    fig.text(0.015, 0, "")
    ax.text(-0.66, yg.mean(), gl[g], transform=ax.get_yaxis_transform(), fontsize=7, color=MUTED, ha="center",
            va="center", linespacing=1.0)
h = [Line2D([], [], color=col if f else "white", marker=mk, mec=col, lw=0, ms=5.5, label=lab)
     for _, lab, col, mk, f in CLS]
fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.02, 0.425), fontsize=7, ncol=3, handletextpad=0.2,
           columnspacing=1.0, borderaxespad=0)
pct = piv.loc[("LAMTOR1 S98", "cat_mTORi")].pct_vs_sites
ax.annotate(f"{pct:.2f}th percentile\nof all sites", (piv.loc[("LAMTOR1 S98", "cat_mTORi")].med_z, y[0]),
            xytext=(-9.3, y[0] + 0.25), fontsize=6.8, color=INK, va="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))

# ---------- (b) turnover reproducibility ----------
ax = aB
ygrid(ax)
FAM = [("ALL", "all"), ("ERK", "ERK"), ("MTOR", "mTOR"), ("AKT", "AKT"), ("CDK1/2/5/9", "CDK")]
ls = h1[h1.analysis == "lab_split"].set_index(["set", "group"])
nml = nm[(nm.set == "random_noise_matched") & (nm.level == "lab")].set_index("group")
SER = [("substrate", "kinase substrates", GREEN, -0.24), ("random", "random sites", GREY, 0.0),
       ("noise", "random sites matched for noise", PURPLE, 0.24)]
for i, (g, gl_) in enumerate(FAM):
    for key, lab, col, off in SER:
        r = nml.loc[g] if key == "noise" else ls.loc[(key, g)]
        ax.bar(i + off, r.rho, width=0.22, color=col, edgecolor="none", zorder=2)
        if not np.isnan(r.lo):
            ax.plot([i + off, i + off], [r.lo, r.hi], color=INK, lw=0.8, zorder=3)
ax.axhline(0, color="#999999", lw=0.8)
ax.set_xticks(range(len(FAM)))
ax.set_xticklabels([f"{l}\n{int(ls.loc[('substrate', g)].n_sites):,}" for g, l in FAM], fontsize=7.8)
ax.text(-0.02, -0.105, "n substrates", transform=ax.transAxes, fontsize=6.3, color=MUTED, ha="right", va="top")
ax.set_ylim(-0.07, 0.6)
ax.set_ylabel("split-by-lab reproducibility\nof response speed (Spearman ρ)")
h = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c, _ in SER]
fig.legend(h, [l for _, l, _, _ in SER], loc="upper left", bbox_to_anchor=(0.66, 0.425), fontsize=7,
           handlelength=1, handletextpad=0.4, borderaxespad=0)

# ---------- (c) ATP consistency ----------
ax = aC
xgrid(ax)
d = atp[atp.shift_tsub_davis.notna()].copy()
d["atp"] = d.atp_site & (d.Km_ATP_uM.notna()) & (d.CP_pred_1mM < 0)
d = d.sort_values("shift_tsub_davis", ascending=False)
yc = np.arange(len(d))
for yi, (_, r) in zip(yc, d.iterrows()):
    if r.atp:
        ax.plot([r.CP_pred_5mM, r.CP_pred_1mM], [yi, yi], color="#d0d0d0", lw=7, solid_capstyle="butt", zorder=1)
        inband = r.CP_pred_5mM <= r.shift_tsub_davis <= r.CP_pred_1mM
        ax.scatter(r.shift_tsub_davis, yi, s=40, color=GREEN if inband else ORANGE, edgecolor="white", lw=0.6,
                   zorder=3)
        ax.text(0.98, yi, f"{imp.loc[r.agent, 'implied_ATP_mM']:.1f} mM", fontsize=7, va="center", ha="left",
                color=INK)
    else:
        ax.scatter(r.shift_tsub_davis, yi, s=40, color=ORANGE, marker="s", edgecolor="white", lw=0.6, zorder=3)
        ax.text(0.98, yi, "n/a", fontsize=7, va="center", ha="left", color=MUTED)
lab = {a: f"{a} → {t}" for a, t, at in zip(d.agent, d.target, d.atp)}
ax.set_yticks(yc)
ax.set_yticklabels([lab[a] for a in d.agent], fontsize=7.6)
ax.axvline(0, color="#999999", lw=0.8)
ax.set_xlim(-3.3, 0.95)
ax.set_ylim(-0.6, len(d) + 0.3)
ax.set_xlabel("potency shift in cells: pEC50 of downstream sites − in-vitro pKd (log10)")
ax.text(0.98, len(d) - 0.45, "implied\n[ATP]", fontsize=6.6, color=MUTED, va="bottom", ha="left")
h = [Line2D([], [], color="#d0d0d0", lw=6, label="predicted for 1–5 mM ATP"),
     Line2D([], [], color=GREEN, marker="o", lw=0, ms=6, label="ATP-site drug (observed)"),
     Line2D([], [], color=ORANGE, marker="s", lw=0, ms=6, label="allosteric MEK inhibitor")]
ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.2, 0.85), fontsize=7, handletextpad=0.4)
med = imp.implied_ATP_mM.median()
fig.text(0.81, 0.13, f"median implied [ATP]:\n{med:.1f} mM (range {imp.implied_ATP_mM.min():.1f}–"
         f"{imp.implied_ATP_mM.max():.1f})\nn = {len(imp)} ATP-site drugs", fontsize=7, color=INK, va="top")

ptitle(fig, aA, "a", "Candidate mTORC1 sites at the lysosome", dx=0.26, dy=0.03)
ptitle(fig, aB, "b", "Response speed reproduces", dx=0.105, dy=0.03)
ptitle(fig, aC, "c", "Cellular potency fits ~2 mM ATP", dx=0.225, dy=0.03)
save(fig, "fig7_caveated")
