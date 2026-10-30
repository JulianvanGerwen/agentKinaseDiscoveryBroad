"""I10 figures. fig1 LOSO AUROC (all/down/up); fig2 CPTAC validation footprint vs kinase vs matched-null; fig3 module/footprint tumour coherence;
fig4 reproducible footprint-feature associations (meta z)."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
RT = "/workspace/results/b_I10/tables/"; FG = "/workspace/results/b_I10/figures/"; CA = "/workspace/cache/b_I10/"
C1, C2, C3, MUTED, INK = '#2a78d6', '#eb6834', '#1baf7a', '#898781', '#0b0b0b'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#c3c2b7', 'axes.labelcolor': INK, 'xtick.color': '#52514e', 'ytick.color': '#52514e',
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb'})
# fig1
A = pd.read_csv(RT + 'T2_loso_auroc.tsv', sep='\t').set_index('footprint').sort_values('aurocOA_all')
fig, ax = plt.subplots(figsize=(6.5, 8)); y = np.arange(len(A))
for off, col, c, lab in [(-0.25, 'aurocOA_all', C1, 'all sites'), (0, 'aurocOA_down', C2, 'down-sites only'), (0.25, 'aurocOA_up', C3, 'up-sites only')]:
    ax.scatter(A[col], y + off, s=22, color=c, label=lab, edgecolor='#fcfcfb', linewidth=1, zorder=3)
ax.axvline(0.5, color=MUTED, lw=1, ls='--'); ax.set_yticks(y); ax.set_yticklabels([f'{i} (n={int(n)})' for i, n in zip(A.index, A.n_pos)])
ax.set_xlabel('Leave-one-study-out AUROC (held-out class vs other-axis experiments)'); ax.legend(frameon=False, loc='lower right')
ax.set_title('Footprints recover held-out studies; inhibitor down-sites transfer better than up-sites', loc='left', fontsize=9)
ax.grid(axis='x', color='#e1e0d9', lw=0.6); fig.tight_layout(); fig.savefig(FG + 'fig1_loso_auroc.png', dpi=160); plt.close()
# fig2
R = pd.read_csv(RT + 'T4_validation.tsv', sep='\t'); R = R[R.variant == 'tc_padj']
tests = ['RAS/RAF hotspot', 'PIK3CA hotspot/PTEN', 'EGFR amp (GBM)', 'TNF/NFkB hallmark', 'PROGENy MAPK (RNA)', 'Immune (ESTIMATE)']
fig, axs = plt.subplots(2, 3, figsize=(12, 7))
for ax, t in zip(axs.flat, tests):
    d = R[R.test == t].sort_values('pooled')
    cols = [C2 if f.startswith('kin:') else (MUTED if f.startswith('ann:') else C1) for f in d.feature]
    ax.barh(range(len(d)), d.pooled - (0.5 if 'AUROC' or True else 0) * (d.test.iloc[0] in ['RAS/RAF hotspot', 'PIK3CA hotspot/PTEN', 'EGFR amp (GBM)']),
            left=0.5 if t in ['RAS/RAF hotspot', 'PIK3CA hotspot/PTEN', 'EGFR amp (GBM)'] else 0, color=cols, height=0.7)
    ax.set_yticks(range(len(d))); ax.set_yticklabels([f.replace('ann:PTM_SEA_PERT-PSP_EGF', 'PTM-SEA EGF (PTMsigDB)') + ('' if pd.isna(p) else f'  p={p:.3f}') for f, p in zip(d.feature, d.p_null)], fontsize=7)
    ax.axvline(0.5 if t in ['RAS/RAF hotspot', 'PIK3CA hotspot/PTEN', 'EGFR amp (GBM)'] else 0, color=MUTED, lw=1)
    ax.set_title(t + (' (pooled AUROC)' if t in ['RAS/RAF hotspot', 'PIK3CA hotspot/PTEN', 'EGFR amp (GBM)'] else ' (mean Spearman)'), loc='left', fontsize=8.5)
fig.text(0.01, 0.005, 'blue = footprint pathway activity (protein-adjusted, p = matched random-site-set null); orange = core-v3 kinase activity; grey = PTMsigDB PTM-SEA EGF', fontsize=7.5, color='#52514e')
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(FG + 'fig2_cptac_validation.png', dpi=160); plt.close()
# fig3
if os.path.exists(RT + 'T8_coherence.tsv'):
    T = pd.read_csv(RT + 'T8_coherence.tsv', sep='\t'); T['kind'] = T.set.str.split(':').str[0]
    S = T.groupby('set').agg(kind=('kind', 'first'), med_z=('z_matched', 'median'), nsig=('p_matched', lambda p: (p < 0.05).sum()))
    fig, ax = plt.subplots(figsize=(6, 3.6)); kinds = ['PSP', 'FP', 'ICA', 'ICAsc']; labs = ['PSP kinase\nsubstrate sets', 'Perturbation\nfootprints', 'ICA modules', 'ICA modules\n(study-centred)']
    for i, k in enumerate(kinds):
        v = S[S.kind == k].nsig.values; ax.scatter(i + np.random.default_rng(0).uniform(-0.2, 0.2, len(v)), v, s=14, color=[C2, C1, C3, C3][i], alpha=0.8, edgecolor='#fcfcfb', linewidth=0.8)
        ax.text(i, 10.6, f'{(v >= 3).mean():.0%} in ≥3', ha='center', fontsize=7.5, color='#52514e')
    ax.set_xticks(range(4)); ax.set_xticklabels(labs); ax.set_ylabel('cohorts (of 10) with signed coherence\n> matched random sets (p<0.05)'); ax.set_ylim(-0.5, 11.2)
    ax.set_title('Do compendium-derived site sets co-vary in tumours?', loc='left', fontsize=9); fig.tight_layout(); fig.savefig(FG + 'fig3_coherence.png', dpi=160); plt.close()
# fig4
if os.path.exists(RT + 'T9_assoc_footprint.tsv'):
    M = pd.read_csv(RT + 'T9_assoc_footprint.tsv', sep='\t'); r = M[M.reproducible]
    top_f = r.groupby('feature').meta_z.apply(lambda z: z.abs().max()).sort_values(ascending=False).index[:15]
    top_s = r.groupby('score').meta_z.apply(lambda z: z.abs().max()).sort_values(ascending=False).index[:20]
    H = M[M.feature.isin(top_f) & M.score.isin(top_s)].pivot(index='score', columns='feature', values='meta_z').reindex(index=top_s, columns=top_f)
    fig, ax = plt.subplots(figsize=(8, 6.5)); v = np.nanmax(np.abs(H.values))
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#f0efec', '#d95926'])
    im = ax.imshow(H.values, cmap=cm, vmin=-v, vmax=v, aspect='auto')
    rep = M.set_index(['score', 'feature']).reproducible
    for i, s in enumerate(H.index):
        for j, f in enumerate(H.columns):
            if rep.get((s, f), False): ax.text(j, i, '•', ha='center', va='center', fontsize=8, color=INK)
    ax.set_xticks(range(len(H.columns))); ax.set_xticklabels(H.columns, rotation=60, ha='right', fontsize=7); ax.set_yticks(range(len(H.index))); ax.set_yticklabels(H.index, fontsize=7)
    plt.colorbar(im, ax=ax, label='meta z (purity-adjusted, 10 cohorts)', shrink=0.6)
    ax.set_title('Tumour perturbation-mimicry vs genomic/clinical features (• reproducible: q<0.05, ≥3 cohorts p<0.05)', loc='left', fontsize=8.5)
    fig.tight_layout(); fig.savefig(FG + 'fig4_mimicry_assoc.png', dpi=160); plt.close()
print('figures done')
