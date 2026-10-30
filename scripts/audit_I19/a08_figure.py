"""Key figure: cross-study replication (log odds ratio vs per-direction base rate) by |log2FC| bin, up vs down,
kinase inhibitors (left) vs activating stimulations (right). Error bars: cluster bootstrap over studies (95%)."""
import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = '/workspace/results/audit_I19/'
Ri = pd.read_csv(OUT + 'cache/pairs_raw_disjoint.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
Rs = pd.read_csv(OUT + 'cache/pairs_stim.tsv.gz', sep='\t', dtype={'si': str, 'sj': str})
Rs = Rs[Rs.group.isin(['GF_RTK', 'PPase_inh', 'DNA_damage', 'GPCR', 'TNF', 'PMA'])]
rng = np.random.default_rng(0)
def est(R, d, b, nb=500):
    sp = R[(R.dir == d) & (R.bin == b)].groupby(['si', 'sj']).lor.mean()
    si, sj = sp.index.get_level_values(0).values, sp.index.get_level_values(1).values
    st = np.unique(np.r_[si, sj]); bs = []
    for _ in range(nb):
        w = pd.Series(rng.choice(st, len(st))).value_counts(); ww = w.reindex(si).fillna(0).values * w.reindex(sj).fillna(0).values
        if ww.sum(): bs.append(np.average(sp.values, weights=ww))
    return sp.mean(), *np.percentile(bs, [2.5, 97.5])
bins = ['0.5-1', '1-2', '2-99']; lab = ['0.5–1', '1–2', '>2']
col = {'down': '#2a78d6', 'up': '#eb6834'}; mk = {'down': 'o', 'up': 's'}
fig, axs = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
rows = []
for ax, (title, R) in zip(axs, [('Kinase inhibitors (23 classes with ≥2 studies)', Ri), ('Activating stimulations (GF/RTK, TNF, DNA damage,\nGPCR, PPase-inh, PMA)', Rs)]):
    for k, d in enumerate(['down', 'up']):
        e = np.array([est(R, d, b) for b in bins]); x = np.arange(3) + (k - .5) * 0.08
        ax.errorbar(x, e[:, 0], yerr=[e[:, 0] - e[:, 1], e[:, 2] - e[:, 0]], color=col[d], marker=mk[d], ms=8, lw=2, capsize=0, elinewidth=1.2)
        ax.text(2.12, e[2, 0], 'decreases' if d == 'down' else 'increases', color='#0b0b0b', va='center', fontsize=9)
        for b, r in zip(bins, e): rows.append(dict(panel=title.split(' (')[0], dir=d, bin=b, lor=r[0], lo=r[1], hi=r[2]))
    ax.axhline(0, color='#999', lw=.8); ax.set_xticks(range(3), lab); ax.set_xlim(-.3, 2.7)
    ax.set_xlabel('|log2FC| in discovery experiment'); ax.set_title(title, fontsize=10, loc='left')
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    ax.grid(axis='y', color='#eee')
axs[0].set_ylabel('cross-study replication\n(log odds ratio vs base rate)')
fig.tight_layout(); fig.savefig(OUT + 'fig_mirror_inhibitor_vs_stimulation.png', dpi=160)
pd.DataFrame(rows).to_csv(OUT + 'T_fig_mirror.tsv', sep='\t', index=False); print(pd.DataFrame(rows).round(2))
