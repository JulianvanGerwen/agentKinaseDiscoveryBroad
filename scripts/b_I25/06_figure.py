import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
T = '/workspace/results/b_I25/tables/'; REV, PRI, GRY, INK = '#2a78d6', '#eb6834', '#8a8984', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})
fig, ax = plt.subplots(1, 3, figsize=(15, 5.4), gridspec_kw=dict(width_ratios=[1.1, 1.1, 1]))
# A
H = pd.read_csv(T + 'T1_H1_enrichment.tsv', sep='\t'); H = H[H.kind == 'inh']
sets = ['CC_Mup', 'CC_Mdown', 'CC_MKI', 'CC_union', 'CC_inter', 'S_MET24', 'S_AMPKmotif', 'S_CULT3D', 'S_HIPPO', 'S_SERUM', 'S_OIS']
H = H.set_index('set').loc[sets]; y = np.arange(len(sets))
for off, c, k, lab in [(-.15, REV, 'rev', 'reverse (up) responders'), (.15, PRI, 'prim', 'primary (down) responders')]:
    v = H[f'logOE_{k}']; ax[0].errorbar(v, y + off, xerr=[v - H[f'{k}_lo'], H[f'{k}_hi'] - v], fmt='o', ms=5, color=c, lw=1.5, label=lab)
ax[0].axvline(0, color=GRY, lw=1); ax[0].set_yticks(y, sets); ax[0].invert_yaxis()
ax[0].set_xlabel('log(observed / matched expected) set membership'); ax[0].set_title('A  H1/H4: axis-set enrichment (inhibitors)', loc='left', fontsize=10)
ax[0].legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
# B
R = pd.read_csv(T + 'T2_H2_removal.tsv', sep='\t'); R = R[(R.scope == 'class_xstudy_all') & (R.bin == '1-99')].set_index('removed')
order = ['none'] + [s for s in sets if s in R.index]; R = R.loc[order]; y = np.arange(len(order))
ax[1].scatter(R.down, y, color=PRI, s=30, label='primary (down) LOR', zorder=3); ax[1].scatter(R.up, y, color=REV, s=30, label='reverse (up) LOR', zorder=3)
rnd = R.up.iloc[0] + R.rnd_delta_up_mean.fillna(0)
ax[1].scatter(rnd.iloc[1:], y[1:], marker='|', s=120, color=GRY, label='reverse after matched random removal', zorder=2)
for yy, (u, d) in enumerate(zip(R.up, R.down)): ax[1].plot([u, d], [yy, yy], color='#d0cfca', lw=1, zorder=1)
ax[1].set_yticks(y, ['all sites'] + [f'minus {s}' for s in order[1:]]); ax[1].invert_yaxis()
ax[1].set_xlabel('cross-study replication LOR (same class, |log2FC|>1)'); ax[1].set_title('B  H2: removing axis sites does not close the gap', loc='left', fontsize=10)
ax[1].legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.4, -0.13), ncol=2)
# C
Q = pd.read_csv(T + 'T4_H4_discordance_R2.tsv', sep='\t'); Q = Q[(Q.outcome == 'up') & Q.adjusted_for.isna()]
Q = Q[~Q.predictor.str.contains('down_LOR')].sort_values('R2'); y = np.arange(len(Q))
ax[2].barh(y, Q.R2, color=[REV if p.startswith('JOINT') else '#9fc2ea' for p in Q.predictor], height=.6)
ax[2].errorbar(Q.R2, y, xerr=[(Q.R2 - Q.lo).clip(lower=0), (Q.hi - Q.R2).clip(lower=0)], fmt='none', color=INK, lw=1)
ax[2].set_yticks(y, Q.predictor.str.replace('d_', '|Δ| ')); ax[2].axvline(.10, color=GRY, ls='--', lw=1)
ax[2].set_xlabel('R² of per-pair reverse replication (LOR_up)'); ax[2].set_title('C  H4: variance in cross-lab reverse replication', loc='left', fontsize=10)
fig.tight_layout(); fig.savefig('/workspace/results/b_I25/figures/fig1_reverse_layer_state_axes.png', dpi=160)
