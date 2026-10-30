"""fig3: Finding K1 - KPNA4 S60 is a CK2 site. Re-plots saved tables only."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from style import *

# ---------- data ----------
rank = pd.read_csv(f"{R}/audit_D/tables/kpna4_ck2i_rank.tsv", sep="\t")
cx = rank[(rank.agent == "CX4945") & rank.z.notna()].copy()
mapk = pd.read_csv(f"{R}/val_I03/tables/t3_mapk_concordance.tsv", sep="\t")
erk = mapk[(mapk.pathway_down == True) & mapk.KPNA4_S60.notna()].copy()
assert len(erk) == 25
men = pd.read_csv(f"{R}/audit_F/task2_menyhart2023_module_sites.tsv", sep="\t")
men = men[~men.log2_MoverL_inhVsDMSO_WT.isin(["NOT_IN_TABLE", "NOT_QUANTIFIED"])].copy()
men = men[men.multiplicity.astype(float) == 1]
for c in ["log2_MoverL_inhVsDMSO_WT", "log2_HoverL_TMinh_vs_DMSO", "log2_HoverM_rescue", "FDR_H_M", "FDR_M_L"]:
    men[c] = pd.to_numeric(men[c], errors="coerce")
psp = pd.read_csv("/workspace/docs/briefing/figures/src_data/derived_menyhart2023_psp_ck2_sites.tsv", sep="\t")

fig = plt.figure(figsize=(7.2, 5.4))
axA = fig.add_axes([0.11, 0.39, 0.33, 0.52])
axB = fig.add_axes([0.6, 0.39, 0.385, 0.52])
axC = fig.add_axes([0.03, 0.035, 0.94, 0.15])

# ---------- (a) compendium ----------
ax = axA
ygrid(ax)
ax.axhline(0, color="#999999", lw=0.8)
ax.axhline(-2, color=GREY, lw=0.8, ls=(0, (3, 2)))
ax.text(1.62, -1.85, "z = −2", fontsize=7.5, color=MUTED, ha="right", va="bottom")
# CX-4945 per cell line
lab = {"HELA": "HeLa (mitotic)", "NTERA2": "NTERA2", "HL60": "HL60"}
for _, r in cx.iterrows():
    ax.scatter(0, r.z, s=46, color=GREEN, edgecolor="white", lw=0.8, zorder=3)
    pct = r.pct
    ax.text(0.12, r.z, f"{lab[r.cell]}  z {r.z:.1f}".replace("-", "−") + f"\nrank {int(r['rank'])} of {int(r.n_sites):,} sites ({pct:.1f}th pct)",
            fontsize=7, va="center", ha="left", color=INK, linespacing=1.05)
# ERK-pathway verified units
rng = np.random.default_rng(1)
x = 1.25 + rng.uniform(-0.17, 0.17, len(erk))
ax.scatter(x, erk.KPNA4_S60, s=16, color=GREY, edgecolor="white", lw=0.5, zorder=3)
med = erk.KPNA4_S60.median()
ax.plot([1.03, 1.47], [med, med], color=INK, lw=1.4, zorder=4)
ax.text(1.25, 2.75, f"median {med:+.2f}\n{(erk.KPNA4_S60 < -2).sum()} of {len(erk)} below −2",
        fontsize=7.5, ha="center", va="bottom", color=INK)
ax.set_xlim(-0.3, 1.65)
ax.set_ylim(-7.8, 4.2)
ax.set_xticks([0, 1.25])
ax.set_xticklabels(["CK2 inhibitor\nCX-4945\n(3 cell lines)", "ERK pathway\ninhibited\n(25 experiments)"])
ax.set_ylabel("change in KPNA4 S60 phosphorylation\n(z-score within experiment)")

# ---------- (b) Menyhart 2023 chemical-genetic rescue ----------
ax = axB
ygrid(ax)
ax.axhline(0, color="#999999", lw=0.8)
arms = [("SGC-CK2-1", "4h"), ("SGC-CK2-1", "24h"), ("CX-4945", "4h"), ("CX-4945", "24h")]
dx = 0.2
for i, (inh, t) in enumerate(arms):
    xw, xr = i - dx, i + dx
    # known CK2 substrates (PSP) median
    p = psp[(psp.inhibitor == inh) & (psp.time == t)]
    pm, ph = p["log2_M_L_WT_inh_vs_DMSO"].median(), p["log2_H_L_TM_inh_vs_DMSO"].median()
    ax.plot([xw, xr], [pm, ph], color=GREY, lw=1.3, zorder=2)
    ax.scatter([xw, xr], [pm, ph], marker="D", s=22, color=GREY, edgecolor="white", lw=0.6, zorder=3)
    # fast CK2 set members other than KPNA4
    m = men[(men.inhibitor == inh) & (men.time == t) & (men.gene != "KPNA4")]
    for _, r in m.iterrows():
        if np.isnan(r.log2_HoverL_TMinh_vs_DMSO):
            ax.scatter(xw, r.log2_MoverL_inhVsDMSO_WT, s=14, facecolor="white", edgecolor=GREEN, lw=0.9, zorder=3)
            continue
        ax.plot([xw, xr], [r.log2_MoverL_inhVsDMSO_WT, r.log2_HoverL_TMinh_vs_DMSO], color=GREEN, lw=0.8, alpha=0.7, zorder=2)
        ax.scatter([xw, xr], [r.log2_MoverL_inhVsDMSO_WT, r.log2_HoverL_TMinh_vs_DMSO], s=14, color=GREEN,
                   edgecolor="white", lw=0.5, zorder=3)
    # KPNA4 S60
    k = men[(men.inhibitor == inh) & (men.time == t) & (men.gene == "KPNA4")].iloc[0]
    ax.plot([xw, xr], [k.log2_MoverL_inhVsDMSO_WT, k.log2_HoverL_TMinh_vs_DMSO], color=PINK, lw=2.0, zorder=4)
    ax.scatter([xw, xr], [k.log2_MoverL_inhVsDMSO_WT, k.log2_HoverL_TMinh_vs_DMSO], s=50, color=PINK,
               edgecolor="white", lw=0.8, zorder=5)
    fdr = k.FDR_H_M
    ax.text(i, -4.55, f"rescue\n{k.log2_HoverM_rescue:+.2f}\nFDR {fdr:.0e}".replace("e-0", "e-"),
            ha="center", va="top", fontsize=6.8, color=INK, linespacing=1.0)
ax.axvline(1.5, color=GRID, lw=1)
ax.set_xticks(range(4))
ax.set_xticklabels([f"{a}\n{t.replace('h', ' h')}" for a, t in arms], fontsize=7.8)
ax.tick_params(axis="x", which="major", length=0, pad=13)
ax.set_xticks([i + s_ for i in range(4) for s_ in (-dx, dx)], minor=True)
ax.set_xticklabels(["WT", "res."] * 4, minor=True, fontsize=6.8, color=MUTED)
ax.tick_params(axis="x", which="minor", length=2, pad=2)
ax.set_xlim(-0.55, 3.55)
ax.set_ylim(-6.1, 2.3)
ax.set_yticks([-4, -3, -2, -1, 0, 1, 2])
ax.set_ylabel("change in phosphorylation\n(log2 fold change vs DMSO)")
h = [Line2D([], [], color=PINK, marker="o", lw=2, ms=6.5, label="KPNA4 S60"),
     Line2D([], [], color=GREEN, marker="o", lw=0.8, ms=4, label="other fast-CK2-set sites (5)"),
     Line2D([], [], color=GREY, marker="D", lw=1.3, ms=4, label="known CK2 substrates (median)")]
ax.legend(handles=h, loc="upper left", bbox_to_anchor=(-0.22, -0.2), ncol=2, fontsize=7, handlelength=1.6,
          columnspacing=0.8, borderaxespad=0)
fig.text(0.515, 0.215, "WT = wild-type CK2α;  res. = inhibitor-resistant CK2α (both + inhibitor)",
         fontsize=7, color=MUTED, ha="left", va="top")

# ---------- (c) sequence ----------
ax = axC
ax.axis("off")
seq = "HEDICEDSDIDGDYR"
c0 = 7  # S60
cw = 0.032
x0 = 0.5 - cw * len(seq) / 2 - 0.16
for j, aa in enumerate(seq):
    pos = j - c0
    xx = x0 + j * cw
    if pos == 0:
        fc, tc = PINK, "white"
    elif aa in "DE":
        fc, tc = GREEN, "white"
    else:
        fc, tc = "#f2f2f2", INK
    ax.add_patch(plt.Rectangle((xx, 0.35), cw * 0.92, 0.42, transform=ax.transAxes, facecolor=fc,
                               edgecolor="none"))
    ax.text(xx + cw * 0.46, 0.56, aa, transform=ax.transAxes, ha="center", va="center", fontsize=10,
            fontweight="bold", color=tc, family="DejaVu Sans Mono")
    ax.text(xx + cw * 0.46, 0.18, ("S60" if pos == 0 else f"{pos:+d}"), transform=ax.transAxes, ha="center",
            va="center", fontsize=6.5, color=MUTED)
# bracket S..+3
xs = x0 + c0 * cw + cw * 0.46
xe = x0 + (c0 + 3) * cw + cw * 0.46
ax.plot([xs, xs, xe, xe], [0.82, 0.92, 0.92, 0.82], transform=ax.transAxes, color=INK, lw=0.9)
ax.text((xs + xe) / 2, 0.97, "CK2 motif: S-x-x-D/E", transform=ax.transAxes, ha="center", va="bottom",
        fontsize=8.5, color=INK)
xl = x0 + len(seq) * cw + 0.03
ax.text(xl, 0.62, "KPNA4 (importin-α3) around S60", transform=ax.transAxes, fontsize=8.5, va="center", color=INK)
ax.add_patch(plt.Rectangle((xl, 0.18), 0.018, 0.2, transform=ax.transAxes, facecolor=GREEN, edgecolor="none"))
ax.text(xl + 0.025, 0.28, "acidic residue (D/E)", transform=ax.transAxes, fontsize=7.5, va="center", color=INK)
ax.add_patch(plt.Rectangle((xl + 0.19, 0.18), 0.018, 0.2, transform=ax.transAxes, facecolor=PINK, edgecolor="none"))
ax.text(xl + 0.215, 0.28, "phospho-site", transform=ax.transAxes, fontsize=7.5, va="center", color=INK)

ptitle(fig, axA, "a", "Compendium: CK2 vs ERK inhibition", dx=0.085)
ptitle(fig, axB, "b", "Chemical-genetic rescue (public data)", dx=0.1)
fig.text(0.01, 0.2, "c", fontsize=11, fontweight="bold", color=INK, va="bottom")
save(fig, "fig3_kpna4_ck2")
