"""Main figure (4 panels) from tables T1b, T3b, T4/T4c, T6."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
T = '/workspace/results/b_I05/tables/'; F = '/workspace/results/b_I05/figures/'
PRI, REV, GRY = '#2a6fdb', '#c96a00', '#7a7a7a'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#999', 'axes.grid': True,
                     'grid.color': '#eee', 'grid.linewidth': 0.6})
fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
# A: reverse-share slope (all responders and adjacent-confirmed)
a = ax[0, 0]
S = pd.read_csv(T + 'T1b_revshare_summary.tsv', sep='\t'); S = S[(S.thr == 1.0) & S.family.isin(['inhibitor', 'stimulus'])]
C = pd.read_csv(T + 'T1c_revshare_confirmed_per_series.tsv', sep='\t', dtype={'study': str})
y = 0; labs = []
for f in ['inhibitor', 'stimulus']:
    r = S[S.family == f].iloc[0]
    a.errorbar(r.slope, y, xerr=[[r.slope - r.lo], [r.hi - r.slope]], fmt='o', color='#333', ms=7, capsize=3, lw=2)
    a.plot([r.null_mean - 2 * r.null_sd, r.null_mean + 2 * r.null_sd], [y - .25] * 2, color=GRY, lw=3, alpha=.6)
    labs.append(f'{f.replace("us","i").replace("inhibitor","inhibitors")}, all responders\n({int(r.n_studies)} studies)'); y += 1
    st = C[C.family == f].groupby('study').slope_conf.mean().dropna()
    a.scatter(st, np.full(len(st), y) + np.random.default_rng(0).uniform(-.12, .12, len(st)), color='#333', s=14, alpha=.7)
    a.plot(st.mean(), y, 'D', color='k', ms=6)
    labs.append(f'{f.replace("us","i").replace("inhibitor","inhibitors")}, confirmed responders\n(>=3 tp; {len(st)} studies; dots=studies)'); y += 1
a.axvline(0, color='#555', lw=1); a.set_yticks(range(y)); a.set_yticklabels(labs); a.invert_yaxis()
a.set_xlabel('slope of logit(reverse share of responders) over time\n(>0: reverse layer grows later)')
a.set_title('A  Reverse-direction responders do not appear later', loc='left', fontweight='bold')
a.text(0.99, 0.02, 'grey bar: shuffled-time null ±2 SD', transform=a.transAxes, ha='right', fontsize=7.5, color=GRY)
# B: portability by onset category
b = ax[0, 1]
P = pd.read_csv(T + 'T3b_portability_persistence_diffcell.tsv', sep='\t'); P = P[(P.kind == 'all') & (P.bin == '1-2+2-99')]
cats = [('early_sus', 'early, sustained'), ('early_trans', 'early, transient'), ('late', 'late onset')]
for i, (c, lab) in enumerate(cats):
    for k, (p, col) in enumerate([('primary', PRI), ('reverse', REV)]):
        r = P[P.contrast == f'{p}|{c}'].iloc[0]; x = i + (k - .5) * .3
        b.errorbar(x, r.est, yerr=[[r.est - r.lo], [r.hi - r.est]], fmt='o', color=col, ms=7, capsize=3, lw=2, label=p if i == 0 else None)
b.set_xticks(range(3)); b.set_xticklabels([l for _, l in cats]); b.axhline(0, color='#555', lw=1)
b.set_ylabel('cross-study, cross-cell replication LOR\n(|peak log2FC|>1; vs per-exp base rate)')
b.set_title('B  Portability tracks persistence, not onset', loc='left', fontweight='bold'); b.legend(frameon=False, loc='upper right')
# C: motif-matched vs unmatched
c_ = ax[1, 0]
M = pd.read_csv(T + 'T4_direct_portability.tsv', sep='\t'); M = M[(M.split == 'motif') & (M.bin == '1-2+2-99') & M.prim.isin(['primary', 'reverse'])]
PL = pd.read_csv(T + 'T4c_placebo_kinase_sets.tsv', sep='\t')
xt = []
for i, kind in enumerate(['stim', 'inh']):
    for k, (p, col) in enumerate([('primary', PRI), ('reverse', REV)]):
        r = M[(M.kind == kind) & (M.prim == p)].iloc[0]; x = i + (k - .5) * .3
        c_.errorbar(x, r.D_est, yerr=[[r.D_est - r.lo], [r.hi - r.D_est]], fmt='o', color=col, ms=7, capsize=3, lw=2, label=p if i == 0 else None)
        pl = PL[(PL.kind == kind) & (PL.prim == p)].D
        c_.scatter(np.full(len(pl), x + .09), pl, marker='_', s=60, color=GRY, label='random kinase sets' if (i == 0 and k == 0) else None)
c_.set_xticks([0, 1]); c_.set_xticklabels(['activating stimuli', 'kinase inhibitors']); c_.axhline(0, color='#555', lw=1)
c_.set_ylabel('Δ replication LOR\n(target-motif-matched − unmatched sites)')
c_.set_title('C  Direct-motif sites are more portable (stimuli only)', loc='left', fontweight='bold'); c_.legend(frameon=False, fontsize=8)
# D: dose
d = ax[1, 1]
D = pd.read_csv(T + 'T6_dose_potency_summary.tsv', sep='\t'); D = D[D.metric == 'd_pEC50']
for i, r in enumerate(D.itertuples()):
    d.errorbar(i, r.est, yerr=[[r.est - r.lo], [r.hi - r.est]], fmt='o', color='k' if r.bin == 'all' else PRI, ms=7, capsize=3, lw=2)
d.set_xticks(range(len(D))); d.set_xticklabels([f'|Δ| {b}' if b != 'all' else 'all bins' for b in D.bin]); d.axhline(0, color='#555', lw=1)
d.set_ylabel('pEC50(decreases) − pEC50(increases)\n(>0: reverse responses right-shifted)')
d.set_title(f'D  decryptM: reverse less potent ({int(D.n_agents.iloc[0])} drugs)', loc='left', fontweight='bold')
fig.tight_layout(); fig.savefig(F + 'fig1_timing_directness_portability.png', dpi=160)
print('ok')
