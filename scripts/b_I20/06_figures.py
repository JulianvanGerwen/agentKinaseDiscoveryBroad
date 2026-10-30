import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from lib import *
BLUE, ORANGE, AQUA, YELLOW, GRAY = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#8a8a85'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#999', 'axes.grid': False})
F = f'{RES}/figures'
# Fig1 first stage
fs = pd.read_parquet(f'{RES}/tables/T1_first_stage_all_genes.parquet')
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
bins = np.linspace(-3, 3, 50)
ax[0].hist(np.log10(fs[~fs.is_kinase].F_cnv.clip(1e-3)), bins=bins, density=True, color=GRAY, alpha=.6, label='non-kinase genes')
ax[0].hist(np.log10(fs[fs.is_kinase].F_cnv.clip(1e-3)), bins=bins, density=True, histtype='step', lw=2, color=BLUE, label='kinases')
ax[0].axvline(1, color='#333', ls='--', lw=1); ax[0].set_xlabel('log10 first-stage F (protein ~ cis-CNV)'); ax[0].set_ylabel('density')
ax[0].legend(frameon=False); ax[0].set_title('Instrument strength, gene x cohort')
k = fs[fs.is_kinase]
ax[1].scatter(k.t_rna_cnv, k.t_cnv, s=6, color=BLUE, alpha=.4, edgecolor='none')
ax[1].axhline(np.sqrt(10), color='#333', ls='--', lw=1)
ax[1].set_xlabel('t (RNA ~ CNV)'); ax[1].set_ylabel('t (protein ~ CNV)'); ax[1].set_title('Kinases: RNA as intermediate')
fig.tight_layout(); fig.savefig(f'{F}/fig1_first_stage.png', dpi=150)
# Fig2 AUROC comparison
e = pd.read_csv(f'{RES}/tables/T2_eval_all.tsv', sep='\t'); e = e[e.is_kinase]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
for a, (lab, nm) in zip(ax, [('psp', 'PSP substrates (n=380 kinase-cohorts)'), ('mot', 'Kinase-Library motif >=90th pct (n=851)')]):
    cols = [(f'auc_{lab}_rf', 'IV RF\ncis-CNV', BLUE), (f'auc_{lab}_nv', 'naive\ncorrelation', ORANGE),
            (f'auc_{lab}_fo', 'IV focal\nCNV', AQUA), (f'auc_{lab}_nb', 'IV adj.\nneighbours', YELLOW),
            (f'perm_{lab}', 'permuted', GRAY), (f'fake_{lab}_mean', 'fake\nnon-kinase', GRAY)]
    data = [e[c].dropna().values for c, _, _ in cols]
    bp = a.boxplot(data, widths=.55, patch_artist=True, showfliers=False, medianprops=dict(color='#222'))
    for p, (_, _, c) in zip(bp['boxes'], cols): p.set_facecolor(c); p.set_alpha(.55); p.set_edgecolor(c)
    for i, d in enumerate(data): a.text(i + 1, a.get_ylim()[1] if False else 0.5, '', ha='center')
    a.axhline(.5, color='#333', lw=.8, ls=':')
    a.set_xticks(range(1, len(cols) + 1)); a.set_xticklabels([f'{l}\nmean {np.mean(d):.3f}' for (_, l, _), d in zip(cols, data)], fontsize=6.5)
    a.set_ylabel('AUROC (signed t, own sites excluded)'); a.set_title(nm)
fig.tight_layout(); fig.savefig(f'{F}/fig2_auroc_iv_vs_naive.png', dpi=150)
# Fig3 EGFR
M = pd.read_parquet(f'{RES}/tables/T3_meta_kinase_site.parquet', columns=['kgene', 'site', 'sgene', 'cis', 'k', 'z', 'n_pos', 'psp_edge', 'p_meta'])
s = M[(M.kgene == 'EGFR') & ~M.cis & (M.k >= 2)].copy(); s['isY'] = s.site.str.contains('_Y')
top = s.nlargest(20, 'z').iloc[::-1]
fig, ax = plt.subplots(figsize=(6, 4.6))
cols = [BLUE if y else GRAY for y in top.isY]
ax.barh(range(len(top)), top.z, color=cols, height=.7)
ax.set_yticks(range(len(top))); ax.set_yticklabels([f"{g} {t.split('_')[1]}{'  [PSP]' if p else ''}  ({n}/{k} coh.)" for g, t, p, n, k in zip(top.sgene.fillna(top.site.str.split('_').str[0]), top.site, top.psp_edge, top.n_pos, top.k)], fontsize=7)
ax.set_xlabel('meta z (EGFR cis-CNV -> protein-adjusted site)')
ax.set_title(f'EGFR top-20: {top.isY.sum()}/20 pY (blue); background {100*s.isY.mean():.1f}% pY', fontsize=9)
fig.tight_layout(); fig.savefig(f'{F}/fig3_EGFR_top_sites.png', dpi=150)
# Fig4 motif null
r = pd.read_csv(f'{RES}/tables/T4b_motif_null.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(6, 3.2)); x = np.arange(len(r))
ax.bar(x - .2, r.own_motif90, .38, color=BLUE, label='own kinase motif >=90')
ax.bar(x + .2, r.null_mean, .38, yerr=r.null_sd, color=GRAY, label='swapped (other-group) motif, mean+-sd')
ax.set_xticks(x); ax.set_xticklabels([f'{a}\n(n={b})' for a, b in zip(r['filter'], r.n)], fontsize=7)
for i, q in enumerate(r.est_FDR_motif): ax.text(i, r.own_motif90[i] + 4, f'FDR~{q:.2f}', ha='center', fontsize=7)
ax.set_ylabel('IV-significant non-PSP sites'); ax.legend(frameon=False, fontsize=7); ax.set_title('Discovery: motif consistency vs null')
fig.tight_layout(); fig.savefig(f'{F}/fig4_discovery_motif_null.png', dpi=150)
print('ok')
