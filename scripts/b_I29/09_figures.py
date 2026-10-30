"""I29 figures: F1 (3 panels: receptor pY / pERK / ETV output by driver group, CPTAC LUAD), F2 forest (robustness + replication)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from lib import *
COL = {'EGFR': '#2a78d6', 'RAS': '#eb6834', 'WT': '#1baf7a'}; INK = '#0b0b0b'; INK2 = '#52514e'; SURF = '#fcfcfb'
LAB = {'EGFR': 'EGFR-mut\n(n=39)', 'RAS': 'RAS-path\n(n=42)', 'WT': 'Driver-neg\n(n=29)'}
A = load('LUAD')
G = pd.read_csv(OUT + 'tables/T1_gate_ischaemia.tsv', sep='\t')
def eff(o, m='FULL: +ISCH_noERK+country+log excision', g='EGFR_vs_ref'):
    r = G[(G.outcome == o) & (G.model == m) & (G.contrast == g)].iloc[0]; return r
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'font.family': 'DejaVu Sans'})
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.9), facecolor=SURF)
panels = [('stage_RTK', 'A  Receptor pY (EGFR/ERBB2/MET)'), ('stage_ERK', 'B  ERK1/2 activation loop pT/pY'), ('ETVc', 'C  MAPK transcriptional output\n(ETV4, ETV5, SPRED1, DUSP6 mRNA)')]
rng = np.random.default_rng(1)
for ax, (o, t) in zip(axs, panels):
    ax.set_facecolor(SURF)
    y = (A[o] - A[o].mean()) / A[o].std()
    for i, g in enumerate(['WT', 'EGFR', 'RAS']):
        v = y[A.group == g].dropna()
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=INK, lw=1.5), boxprops=dict(color=INK2, lw=1), whiskerprops=dict(color=INK2, lw=1), capprops=dict(color=INK2, lw=1))
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=16, color=COL[g], alpha=0.8, edgecolor=SURF, linewidth=0.6, zorder=3)
    ax.axhline(0, color='#d6d5d0', lw=0.8, zorder=0)
    ax.set_xticks(range(3)); ax.set_xticklabels([LAB[g] for g in ['WT', 'EGFR', 'RAS']])
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    ax.set_title(t, loc='left', fontsize=9.5, color=INK)
    r = eff(o); rr = eff(o, g='RAS_vs_ref')
    ax.text(0.02, 0.02, f"EGFR vs driver-neg: {r.beta:+.2f} SD [{r.lo:+.2f}, {r.hi:+.2f}], p={r.p:.1g}\nRAS vs driver-neg: {rr.beta:+.2f} SD, p={rr.p:.2f}", transform=ax.transAxes, fontsize=7.3, color=INK2, va='bottom')
    ax.set_ylim(min(-3.6, y.min() - 1.2), max(3, y.max() + 0.3))
axs[0].set_ylabel('z-score within CPTAC LUAD')
fig.suptitle('EGFR-mutant LUAD: receptor pY ≈/↑, ERK activation loop ↓, ERK transcriptional output ↑   (adjusted for purity, E2F, ischaemia score, country, cold-ischaemia time)', fontsize=8.5, color=INK2, x=0.01, ha='left', y=0.995)
fig.tight_layout(); fig.savefig(OUT + 'figures/F1_driver_groups.png', dpi=200, facecolor=SURF); plt.close(fig)

# ---- F2 forest
T2 = pd.read_csv(OUT + 'tables/T2_robustness_panel.tsv', sep='\t')
C6 = pd.read_csv(OUT + 'tables/T6_replication_chen.tsv', sep='\t'); C8 = pd.read_csv(OUT + 'tables/T8_replication_cptac_conf.tsv', sep='\t'); C9 = pd.read_csv(OUT + 'tables/T9_other_cohorts.tsv', sep='\t')
def g2(o, m): r = G[(G.outcome == o) & (G.model == m) & (G.contrast == 'EGFR_vs_ref')].iloc[0]; return r
def t2(o, m, c='EGFR_vs_ref'): return T2[(T2.outcome == o) & (T2.model == m) & (T2.contrast == c)].iloc[0]
items = [
 ('Discovery: BCM pERK, base', g2('stage_ERK', 'base'), 'pERK'),
 ('  + ischaemia score + country + cold-isch. time', g2('stage_ERK', 'FULL: +ISCH_noERK+country+log excision'), 'pERK'),
 ('  Broad pTpY peptides (MAPK1+MAPK3), protein-adj.', g2('pERK_broad', 'base'), 'pERK'),
 ('  Broad mono-pY peptides', t2('pERK_monoY', 'base (vs WTall)'), 'pERK'),
 ('  clean EGFR (ex19/L858R) vs clean driver-neg', t2('stage_ERK', 'EGFR ex19/L858R no RAS vs WT no fusion/EGFRamp/ERBB2'), 'pERK'),
 ('  never-smokers only', t2('stage_ERK', 'never-smokers only'), 'pERK'),
 ('  within Vietnam only', g2('stage_ERK', 'within Vietnam'), 'pERK'),
 ('Confirmatory CPTAC LUAD: MAPK3 mono-pY (padj)', C8[(C8.outcome == 'MAPK3:y204_padj') & (C8.model == 'base') & (C8.contrast == 'EGFR_vs_WT')].iloc[0], 'pERK'),
 ('Chen 2020 TW: pERK T202+T185 (T/N), vs non-EGFR', C6[(C6.outcome == 'pERK') & (C6.model.str.startswith('ADC EGFR vs all')) & (C6.contrast == 'EGFR_vs_WT')].iloc[0], 'pERK'),
 ('GBM EGFR-mut: Broad pTpY', C9[(C9.cohort == 'GBM') & (C9.outcome == 'pERK_broad') & (C9.contrast == 'EGFRmut_vs_ref')].iloc[0], 'pERK'),
 ('Discovery: ETV composite, base', g2('ETVc', 'base'), 'ETV'),
 ('  + ischaemia score + country + cold-isch. time', g2('ETVc', 'FULL: +ISCH_noERK+country+log excision'), 'ETV'),
 ('  never-smokers only', t2('ETVc', 'never-smokers only'), 'ETV'),
 ('Chen 2020 TW: ETV composite, vs non-EGFR', C6[(C6.outcome == 'ETVc') & (C6.model.str.startswith('ADC EGFR vs all')) & (C6.contrast == 'EGFR_vs_WT')].iloc[0], 'ETV'),
 ('GBM EGFR-mut: ETV composite', C9[(C9.cohort == 'GBM') & (C9.outcome == 'ETVc') & (C9.contrast == 'EGFRmut_vs_ref')].iloc[0], 'ETV'),
]
fig, ax = plt.subplots(figsize=(8.2, 5.6), facecolor=SURF); ax.set_facecolor(SURF)
for i, (lab, r, k) in enumerate(items):
    yv = len(items) - i
    c = COL['EGFR'] if k == 'pERK' else '#6250d6'
    ax.plot([r.lo, r.hi], [yv, yv], color=c, lw=2, solid_capstyle='round'); ax.scatter([r.beta], [yv], s=40, color=c, edgecolor=SURF, linewidth=1.2, zorder=3)
    ax.text(1.95, yv, f"{r.beta:+.2f}  (n {int(r.n_g)} vs {int(r.n_ref)})", va='center', fontsize=7.5, color=INK2)
ax.set_yticks([len(items) - i for i in range(len(items))]); ax.set_yticklabels([x[0] for x in items], fontsize=7.8)
ax.axvline(0, color=INK2, lw=0.8); ax.set_xlim(-2.0, 3.0); ax.set_xticks([-2, -1, 0, 1]); ax.set_xlabel('EGFR-altered vs reference, SD (95% CI, HC3)')
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.set_title("pERK (blue): lower in discovery, not replicated in 2 independent LUAD cohorts\nETV output (violet)", fontsize=8.5, color=INK, loc="left")
fig.tight_layout(); fig.savefig(OUT + 'figures/F2_robustness_replication_forest.png', dpi=200, facecolor=SURF)
print('ok')
