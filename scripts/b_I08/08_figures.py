#!/usr/bin/env python3
"""I08 figures (PNG) in results/b_I08/figures."""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *
FIG = f'{OUT}/figures'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#888', 'axes.titlesize': 10})
COL = {'abundance': '#eb6834', 'abundance+PTM': '#6250d6', 'PTM': '#2a78d6', 'uncoupled': '#b5b3ab'}
C = pd.read_csv(f'{OUT}/tables/kinase_classification_all_sources.tsv', sep='\t')

# F1 protein vs activating-site coupling
fig, axs = plt.subplots(1, 2, figsize=(11, 5))
for ax, src in zip(axs, ['meanz_tc_padj', 'corev3_tc_padj']):
    x = C[(C.source == src) & C.prot_rp.notna() & C.site_rp.notna()]
    for cls, g in x.groupby('cls'):
        col = COL.get(cls, '#b5b3ab')
        ax.errorbar(g.prot_rp, g.site_rp, xerr=[g.prot_rp - g.prot_lo, g.prot_hi - g.prot_rp], yerr=[g.site_rp - g.site_lo, g.site_hi - g.site_rp],
                    fmt='o', ms=5, color=col, ecolor=col, elinewidth=.6, alpha=.85, label=f'{cls} (n={len(g)})', mec='white', mew=.8)
    lab = x[(x.cls != 'uncoupled') | (x.site_rp.abs() > .25) | (x.prot_rp.abs() > .15)]
    for _, r in lab.iterrows(): ax.annotate(r.kinase, (r.prot_rp, r.site_rp), fontsize=7, xytext=(3, 2), textcoords='offset points', color='#333')
    ax.axhline(0, color='#ccc', lw=.8); ax.axvline(0, color='#ccc', lw=.8); ax.plot([-.3, .7], [-.3, .7], ls=':', color='#999', lw=.8)
    ax.set_xlim(-.3, .45); ax.set_ylim(-.45, .8)
    ax.set_xlabel('activity ~ own protein (meta partial r, 95% CI)'); ax.set_ylabel('activity ~ own activating/A-loop site, protein-residualised')
    ax.set_title(f'{src.replace("_tc_padj", "")} activity (10 cohorts, adj. purity/ESTIMATE/prolif.)'); ax.legend(frameon=False, fontsize=7, loc='upper left')
fig.tight_layout(); fig.savefig(f'{FIG}/F1_protein_vs_site_coupling.png', dpi=180); plt.close(fig)

# F2 distribution of meta coupling by feature + nulls
M = pd.read_csv(f'{OUT}/tables/meta_all.tsv', sep='\t')
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, src in zip(axs, ['meanz_tc_padj', 'corev3_tc_padj']):
    S = M[M.source == src]
    groups = [('mRNA', S[S.ftype == 'rna'].r_p), ('protein\n(raw r)', S[S.ftype == 'prot'].r), ('protein\n(partial)', S[S.ftype == 'prot'].r_p),
              ('other kinases\' protein\n(median r)', S[S.ftype == 'prot'].other_median_r),
              ('canonical site\n(raw)', S[S.ftype == 'canon_raw'].r_p), ('canonical site\n(protein-resid.)', S[S.ftype == 'canon_padj'].r_p),
              ('any own site\n(protein-resid.)', S[S.ftype.isin(['site_padj', 'canon_padj'])].r_p)]
    data = [g.dropna().values for _, g in groups]
    bp = ax.boxplot(data, widths=.55, showfliers=False, patch_artist=True)
    for b in bp['boxes']: b.set(facecolor='#dce8f8', edgecolor='#2a78d6')
    for i, d in enumerate(data):
        ax.scatter(np.random.default_rng(i).normal(i + 1, .07, len(d)), d, s=5, color='#2a78d6', alpha=.35, lw=0)
        ax.text(i + 1, -.5, f'n={len(d)}\nmed {np.median(d):.2f}', ha='center', fontsize=6.5, color='#555')
    ax.axhline(0, color='#999', lw=.8); ax.set_xticks(range(1, len(groups) + 1)); ax.set_xticklabels([g for g, _ in groups], fontsize=6.5, rotation=35, ha='right')
    ax.set_ylabel('meta r with substrate-inferred activity'); ax.set_ylim(-.55, .75); ax.set_title(src)
fig.tight_layout(); fig.savefig(f'{FIG}/F2_coupling_distributions.png', dpi=180); plt.close(fig)

# F3 forest: per-cohort partial r for headline features
L = pd.read_parquet(f'{OUT}/cache/coupling_long.parquet')
heads = [('MAPK1', 'P28482_Y187|padj'), ('PTK2', 'Q05397_Y576|padj'), ('MTOR', 'P42345_S2448|padj'), ('PAK2', 'Q13177_S197|padj'),
         ('IKBKB', 'O14920_S675|padj'), ('ABL1', 'P00519_Y185|padj'), ('EGFR', 'prot'), ('CDK1', 'prot'), ('MAPK1', 'prot'), ('PTK2', 'prot')]
fig, axs = plt.subplots(2, 5, figsize=(14, 6), sharex=True)
for ax, (k, f) in zip(axs.flat, heads):
    for j, (src, col) in enumerate([('meanz_tc_padj', '#2a78d6'), ('corev3_tc_padj', '#eb6834')]):
        g = L[(L.source == src) & (L.kinase == k) & (L.feature == f)].sort_values('cohort')
        if g.empty: continue
        y = np.array([cl.COHORTS.index(c) for c in g.cohort]) + j * .3
        se = 1 / np.sqrt(np.maximum(g.n - 6, 1))
        lo, hi = np.tanh(np.arctanh(g.r_partial) - 1.96 * se), np.tanh(np.arctanh(g.r_partial) + 1.96 * se)
        ax.errorbar(g.r_partial, y, xerr=[g.r_partial - lo, hi - g.r_partial], fmt='o', ms=3.5, color=col, elinewidth=.7, label=src.split('_')[0])
        m = meta(g.r_partial.values, g.n.values, 3)
        ax.axvline(m['r'], color=col, lw=1, ls='--')
    ax.set_yticks(np.arange(10) + .15); ax.set_yticklabels(cl.COHORTS, fontsize=7)
    ax.axvline(0, color='#999', lw=.8); ax.set_title(f'{k}: {f.replace("|padj", " (resid.)")}', fontsize=8.5)
axs.flat[0].legend(frameon=False, fontsize=7)
for ax in axs[1]: ax.set_xlabel('partial r with activity')
fig.tight_layout(); fig.savefig(f'{FIG}/F3_forest_headlines.png', dpi=180); plt.close(fig)

# F4 reporters: best reporter site vs canonical site
W = pd.read_csv(f'{OUT}/tables/reporter_sites.tsv', sep='\t')
W['best'] = W[['r_p__meanz', 'r_p__corev3']].mean(1)
fig, ax = plt.subplots(figsize=(7, 6))
W = W.sort_values('best')
cmap = {'known activating/A-loop': '#2a78d6', 'PSP regulatory (other)': '#1baf7a', 'novel (no PSP regulatory annotation)': '#eb6834'}
for i, (_, r) in enumerate(W.iterrows()):
    ax.barh(i, r.best, color=cmap[r.status], height=.7, alpha=.9 if r.both_pass else .45)
    if np.isfinite(r.canon_rp__corev3) and r.site != r.canonical_site:
        ax.plot(r.canon_rp__corev3, i, marker='|', color='k', ms=8)
ax.set_yticks(range(len(W))); ax.set_yticklabels([f'{k} {s.split("_")[1]}{" *" if b else ""}' for k, s, b in zip(W.kinase, W.site, W.both_pass)], fontsize=6.5)
ax.set_xlabel('meta partial r with own substrate-inferred activity (mean of methods)')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=v, label=k) for k, v in cmap.items()] + [plt.Line2D([], [], marker='|', color='k', ls='', label='canonical site (corev3)')], frameon=False, fontsize=7, loc='lower right')
ax.set_title('Activity-reporter phosphosites on kinases (* = passes in both methods; faded = corev3 only)', fontsize=8.5)
fig.tight_layout(); fig.savefig(f'{FIG}/F4_reporter_sites.png', dpi=180); plt.close(fig)
print('ok')

# F5 Q3: LOCO elastic-net held-out Spearman per feature block
E = pd.read_csv(f'{OUT}/tables/drivers_enet_loco.tsv', sep='\t')
blocks = ['kinase_act', 'prot_reg', 'mutation', 'cnv_reg']
names = {'kinase_act': 'other kinases\' activity', 'prot_reg': 'regulator protein', 'mutation': 'driver mutations', 'cnv_reg': 'regulator CNV'}
tg = sorted(E.target.unique())
fig, ax = plt.subplots(figsize=(10, 4))
cols = ['#2a78d6', '#1baf7a', '#eb6834', '#6250d6']
for i, b in enumerate(blocks):
    m = E[E.block == b].groupby('target').rho.median().reindex(tg)
    ax.bar(np.arange(len(tg)) + (i - 1.5) * .2, m.values, width=.18, color=cols[i], label=names[b])
ax.set_xticks(range(len(tg))); ax.set_xticklabels(tg, rotation=40, ha='right', fontsize=8); ax.axhline(0, color='#999', lw=.8)
ax.set_ylabel('held-out-cohort Spearman (median of 10 LOCO folds)'); ax.legend(frameon=False, fontsize=7, ncol=4)
ax.set_title('Q3: what predicts activity of abundance-uncoupled, variable kinases (elastic net, train 9 cohorts / test 1)', fontsize=9)
fig.tight_layout(); fig.savefig(f'{FIG}/F5_drivers_loco_enet.png', dpi=180); plt.close(fig)
print('ok5')
