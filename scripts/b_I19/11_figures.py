"""Figures for I19. Palette: blue=down (#2a78d6), orange=up (#eb6834), gray=null; text in neutral ink."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import *
OUT = "/workspace/results/b_I19/"; F = OUT + 'figures/'
UP, DN, NUL, INK, MUT = '#eb6834', '#2a78d6', '#9a9a94', '#2b2b28', '#6b6b66'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUT, 'axes.labelcolor': INK, 'xtick.color': MUT, 'ytick.color': MUT, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb'})
# Fig 1: replication vs discovery magnitude (pooled), three variants
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
bins = ['0.5-1', '1-2', '2-99']; lab = ['0.5–1', '1–2', '>2']
for ax, (suf, title) in zip(axs, [('', 'within-context centred'), ('_raw', 'uncentred vs DMSO'), ('_lowdown', 'excl. top-25% "massive-down" exps')]):
    R = pd.read_csv(OUT + f'tables/T3_replication_pairs{suf}.tsv', sep='\t')
    for d, c in [('down', DN), ('up', UP)]:
        y, yn = [], []
        for b in bins:
            q = R[(R.dir == d) & (R.bin == b)]
            y.append(np.average(q.rep, weights=q.n_disc)); qq = q.dropna(subset=['rep_null']); yn.append(np.average(qq.rep_null, weights=qq.n_disc))
        ax.plot(lab, y, '-o', color=c, lw=2, ms=6, label=f'{d}: observed')
        ax.plot(lab, yn, '--o', color=c, lw=1.5, ms=5, alpha=0.5, label=f'{d}: matched-random null')
    ax.set_title(title, color=INK, fontsize=9); ax.set_xlabel('|study-level log2FC| in discovery study')
    ax.grid(axis='y', color='#e6e5e0', lw=0.6)
axs[0].set_ylabel('P(same-direction response\nin another study)'); axs[0].legend(frameon=False, fontsize=7)
fig.suptitle('Down-regulated sites replicate across studies in proportion to effect size; up-regulated sites do not', color=INK, fontsize=10)
fig.tight_layout(); fig.savefig(F + 'fig1_up_vs_down_replication.png', dpi=200); plt.close()
# Fig 2: pair-type decomposition
g = pd.read_csv(OUT + 'tables/T3_pairtype_summary.tsv', sep='\t', index_col=0)
order = ['same_class_diff_drug_same_cell', 'same_drug_diff_cell', 'same_class_diff_drug_diff_cell', 'same_class_diff_study', 'NULL_diff_class_same_cell', 'NULL_diff_class_diff_cell']
names = ['same class, other drug,\nsame cell line', 'same drug,\nother cell line', 'same class, other drug,\nother cell line', 'same class,\nother study', 'NULL: other class,\nsame cell line', 'NULL: other class,\nother cell line']
g = g.loc[order]; x = np.arange(len(order)); w = 0.38
fig, ax = plt.subplots(figsize=(8, 3.4))
ax.bar(x - w / 2 - 0.01, g.ex_up, w, color=UP, label='up'); ax.bar(x + w / 2 + 0.01, g.ex_dn, w, color=DN, label='down')
ax.text(0.99, 0.62, 'discovery: >2-fold in exp 1; replication: >1.4-fold same direction in exp 2', transform=ax.transAxes, ha='right', fontsize=7, color=MUT)
for i, (u, d) in enumerate(zip(g.ex_up, g.ex_dn)):
    ax.text(i - w / 2, u + 0.006, f'{u:.2f}', ha='center', fontsize=7, color=INK); ax.text(i + w / 2, d + 0.006, f'{d:.2f}', ha='center', fontsize=7, color=INK)
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=7); ax.set_ylabel('replication above background\n(median over experiment pairs)')
ax.legend(frameon=False, fontsize=8); ax.grid(axis='y', color='#e6e5e0', lw=0.6)
ax.set_title('Up-responses are drug-invariant within a cell line but lost across cell lines/studies', color=INK, fontsize=10)
fig.tight_layout(); fig.savefig(F + 'fig2_pairtype_decomposition.png', dpi=200); plt.close()
# Fig 3: self-hyperphosphorylation profiles
P = pd.read_csv(OUT + 'tables/T1_self_hits_profiles.tsv', sep='\t')
sel = ['P31749_S473', 'P31749_T312', 'O96013_S148', 'P17252_T501', 'Q02750_S222', 'P36507_S226', 'P28482_T185', 'P06493_Y15', 'P15056_S365', 'P06239_T395', 'Q13153_S174']
P = P.set_index('site').loc[[s for s in sel if s in P.site.values]]
fig, ax = plt.subplots(figsize=(8, 3.6)); x = np.arange(len(P))
E = experiments(); X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]
for i, (s, r) in enumerate(P.iterrows()):
    xv = X.loc[s]; own = E[E.bset.apply(lambda b: r.gene in b)]; v = xv[own.exp].dropna()
    ax.scatter(np.full(len(v), i - 0.15) + np.random.default_rng(i).uniform(-0.06, 0.06, len(v)), v, s=14, color=UP, zorder=3, label='drug binds this kinase' if i == 0 else None)
    oth = E[~E.bset.apply(lambda b: r.gene in b)]; vo = xv[oth.exp].dropna()
    ax.boxplot([vo.values], positions=[i + 0.18], widths=0.22, showfliers=False, patch_artist=True, boxprops=dict(facecolor='#e6e5e0', color=NUL), medianprops=dict(color=INK), whiskerprops=dict(color=NUL), capprops=dict(color=NUL))
ax.axhline(0, color=MUT, lw=0.6)
ax.set_xticks(x); ax.set_xticklabels([f"{r.gene} {r.residue_pos}" for s, r in P.iterrows()], rotation=40, ha='right', fontsize=7.5)
ax.set_ylabel('log2FC vs DMSO'); ax.scatter([], [], marker='s', color='#e6e5e0', label='all other drugs (box)'); ax.legend(frameon=False, fontsize=7, loc='upper right')
ax.set_title('Kinase own-site phosphorylation under drugs that bind it (dots) vs other drugs', color=INK, fontsize=10)
fig.tight_layout(); fig.savefig(F + 'fig3_self_hyperphosphorylation.png', dpi=200); plt.close()
# Fig 4: (a) RAF paradox, (b) MET/ALK -> Hippo kinases forest
R = pd.read_csv(OUT + 'tables/T2_raf_paradox.tsv', sep='\t')
A = pd.read_csv(OUT + 'tables/T2_compensation_atlas_all.tsv', sep='\t', dtype={'per_study': str})
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw={'width_ratios': [1, 1.3]})
ax = axs[0]
for j, (b, c) in enumerate([('WT', UP), ('V600', DN)]):
    q = R[R.braf == b]; ax.scatter(np.full(len(q), j) + np.random.default_rng(j).uniform(-0.1, 0.1, len(q)), q.MAPK1, color=c, s=22, zorder=3)
    for _, rr in q.iterrows():
        if rr.exp.startswith('32034161') and rr.exp.split('|')[2] in ('LY3009120', 'sorafenib'): ax.annotate(rr.exp.split('|')[2] + ' (type II, RPE1)', (j + 0.12, rr.MAPK1 + (0.35 if 'sora' in rr.exp else -0.45)), fontsize=6.5, color=MUT)
ax.axhline(0, color=MUT, lw=0.6); ax.set_xticks([0, 1]); ax.set_xticklabels(['BRAF-WT context', 'BRAF-V600 context']); ax.set_ylabel('ERK2 (MAPK1) activity, z_core')
ax.set_title('Positive control: RAF-inhibitor paradox', color=INK, fontsize=9)
ax = axs[1]
rows = A[(A.cls == 'MET_ALK') & (A.method == 'corev3') & A.kinase.isin(['LATS1', 'LATS2', 'STK38', 'STK38L'])]
yy = 0; ticks = []
for _, r in rows.iterrows():
    for kv in str(r.per_study).split(';'):
        st, val = kv.split(':'); ax.scatter(float(val), yy, color=NUL, s=10)
    ax.errorbar(r.est, yy, xerr=1.96 * r.se, fmt='D', color=UP, ms=6, lw=2)
    ticks.append(f'{r.kinase} (k={r.k_studies}, n={r.n_exp})'); yy += 1
ax.axvline(0, color=MUT, lw=0.6); ax.set_yticks(range(len(ticks))); ax.set_yticklabels(ticks, fontsize=8)
ax.set_xlabel('activity change (centred z_core)\ndiamond = random-effects estimate ±95% CI; dots = studies')
ax.set_title('MET/ALK inhibitors → Hippo/NDR kinase activity up', color=INK, fontsize=9)
fig.tight_layout(); fig.savefig(F + 'fig4_raf_paradox_and_hippo.png', dpi=200); plt.close()
print('ok')
