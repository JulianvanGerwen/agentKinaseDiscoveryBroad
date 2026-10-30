import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
T = '/workspace/results/b_I22/tables/'; BL, OR, INK, MU = '#2a78d6', '#eb6834', '#0b0b0b', '#8a8984'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': MU, 'axes.labelcolor': INK})
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), gridspec_kw=dict(width_ratios=[1.1, 1.4, 1]))
# A Hijazi decomposition
H = pd.read_csv(T + 'T2_hijazi_decomposition.tsv', sep='\t').query("defn=='dpct' and cen=='cen'").set_index('dir')
types = ['same drug, diff cell', 'diff class, same cell', 'diff class, diff cell']; x = np.arange(3)
for k, (d, c) in enumerate([('primary', BL), ('reverse', OR)]):
    ax[0].bar(x + (k - .5) * .38, H.loc[d, types], .36, color=c, label=d)
ax[0].set_xticks(x, ['same drug\ndiff cell', 'diff class\nsame cell', 'diff class\ndiff cell']); ax[0].axhline(0, color=MU, lw=.8)
ax[0].set_ylabel('motif-profile r'); ax[0].set_title('A  Within Hijazi (shared DMSO per cell, centred)', loc='left', fontsize=9); ax[0].legend(frameon=False)
# B same-cell excess, multiple tests
T1 = pd.read_csv(T + 'T1_H1_cell_identity.tsv', sep='\t').query("defn=='dpct' and cen=='raw'")
T3 = pd.read_csv(T + 'T3_H1_sitelevel_samecell.tsv', sep='\t')
rows = []
for lab, key in [('motif, x-study, diff class', 'A cross'), ('motif, x-study, inh only', 'B cross'), ('motif, x-study, inh vs stim', 'C cross')]:
    for d in ['primary', 'reverse']:
        r = T1[T1.test.str.startswith(key) & (T1.dir == d)].iloc[0]; rows.append((lab, d, r.delta, r.lo, r.hi))
for d in ['primary', 'reverse']:
    r = T3[(T3.kind == 'inh') & (T3.dir == d)].iloc[0]; rows.append(('site LOR, x-study, diff class', d, r.delta, r.lo, r.hi))
for d in ['primary', 'reverse']:
    r = T1[T1.test.str.startswith('F within') & (T1.dir == d)].iloc[0]; rows.append(('motif, within decryptM (indep. ctrl)', d, r.delta, np.nan, np.nan))
    r = T1[T1.test.str.startswith('G within') & (T1.dir == d)].iloc[0]; rows.append(('motif, within Hijazi (shared ctrl)', d, r.delta, np.nan, np.nan))
R = pd.DataFrame(rows, columns=['lab', 'dir', 'd', 'lo', 'hi']); labs = list(dict.fromkeys(R.lab))
for k, (d, c) in enumerate([('primary', BL), ('reverse', OR)]):
    s = R[R.dir == d].set_index('lab').loc[labs]; y = np.arange(len(labs)) + (k - .5) * .3
    ax[1].errorbar(s.d, y, xerr=[s.d - s.lo, s.hi - s.d], fmt='o', color=c, ms=5, lw=1.2, label=d)
ax[1].set_yticks(np.arange(len(labs)), labs); ax[1].invert_yaxis(); ax[1].axvline(0, color=MU, lw=.8)
ax[1].set_xlabel('same-cell minus diff-cell similarity (95% study bootstrap CI)'); ax[1].set_title('B  Is the reverse layer a property of the cell line?', loc='left', fontsize=9)
# C H2
T5 = pd.read_csv(T + 'T5_H2_H4_baseline_vs_profile.tsv', sep='\t').query("family==False")
T5['lab'] = T5.src.map({'rna': 'RNA', 'prot': 'proteome', 'chronos': 'Chronos'}) + np.where(T5.within, ' (within study)', ' (all units)')
labs = list(dict.fromkeys(T5.lab))
for k, (d, c) in enumerate([('primary', BL), ('reverse', OR)]):
    s = T5[T5.dir == d].set_index('lab').loc[labs]; y = np.arange(len(labs)) + (k - .5) * .3
    ax[2].errorbar(s.null_mean, y, xerr=1.96 * s.null_sd, fmt='|', color=MU, lw=1)
    ax[2].plot(s.mean_rho, y, 'o', color=c, ms=5, label=d)
ax[2].set_yticks(np.arange(len(labs)), labs); ax[2].invert_yaxis(); ax[2].axvline(0, color=MU, lw=.8)
ax[2].set_xlabel('mean Spearman ρ, baseline vs motif profile\n(grey: permuted-cell null ±1.96 sd)'); ax[2].set_title('C  Does baseline predict it?', loc='left', fontsize=9)
plt.tight_layout(); plt.savefig('/workspace/results/b_I22/figures/fig1_reverse_layer_cell_identity.png', dpi=170)
