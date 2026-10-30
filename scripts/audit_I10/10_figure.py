"""audit_I10 figure: (A) median tumour coherence under stress tests vs matched null & PSP; (B) per-cohort r(footprint score, cold-ischaemia score)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
RA = '/workspace/results/audit_I10/'
A = pd.read_csv(RA + 'A1_coherence.tsv', sep='\t'); I = pd.read_csv(RA + 'A5_ischaemia.tsv', sep='\t'); I = I[I.cohort != 'ALL']; I['value'] = I.value.astype(float)
sets = ['FP:KI_MEK', 'FP:EGF', 'FP:TNF', 'FP:DNA_damage', 'FP:phosphatase_inh', 'PSP:MAPK1', 'PSP:CDK1', 'PSP:AKT1']
lab = ['MEKi', 'EGF', 'TNF', 'DNA dmg', 'pervanadate', 'PSP MAPK1', 'PSP CDK1', 'PSP AKT1']
med = lambda c: A[A.cond == c].groupby('set').coherence.median().reindex(sets)
isr = I[I.what == 'coherence after regressing ischaemia score'].groupby('set').value.median()
cols = {'base (protein-adjusted)': (med('base'), '#2a78d6'), '+ co-abundant/STRING pairs removed + top-5 PCs regressed': (med('pc5+coab+string'), '#eb6834'),
        '+ cold-ischaemia score regressed': (pd.Series({s: isr.get(s.split(':')[1], np.nan) if s.startswith('FP') else np.nan for s in sets}), '#1baf7a')}
nullB = A[A.cond == 'base'].groupby('set').nullB_mean.median().reindex(sets)
fig, ax = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={'width_ratios': [1.5, 1]}, facecolor='#fcfcfb')
x = np.arange(len(sets)); w = 0.26
for k, (nm, (v, c)) in enumerate(cols.items()):
    ax[0].bar(x + (k - 1) * w, v.values, w - 0.03, color=c, label=nm)
ax[0].scatter(x, nullB.values, marker='_', s=300, color='#0b0b0b', zorder=5, label='matched random-site null (cov×SD×communality×residue×P+1)')
ax[0].set_xticks(x); ax[0].set_xticklabels(lab, rotation=30, ha='right'); ax[0].set_ylabel('signed mean cross-protein pairwise r\n(median of 10 CPTAC cohorts)')
ax[0].legend(frameon=False, fontsize=8, loc='upper right'); ax[0].set_title('A  Footprint coherence survives confound removal but shrinks ~2x', loc='left', fontsize=10)
fp = ['KI_MEK', 'EGF', 'TNF', 'phosphatase_inh', 'GPCR_ligand', 'DNA_damage', 'metformin']
r = I[I.what == 'r(score, ischaemia score)']
for j, g in enumerate(fp):
    v = r[r.set == g].value.values
    ax[1].scatter(np.full(len(v), j) + np.random.default_rng(j).uniform(-.12, .12, len(v)), v, s=28, color='#2a78d6', edgecolor='#fcfcfb', linewidth=1)
    ax[1].plot([j - .25, j + .25], [np.median(v)] * 2, color='#0b0b0b', lw=2)
ax[1].axhline(0, color='#52514e', lw=.8); ax[1].set_xticks(range(len(fp))); ax[1].set_xticklabels(['MEKi', 'EGF', 'TNF', 'pervanadate', 'GPCR', 'DNA dmg', 'metformin'], rotation=30, ha='right')
ax[1].set_ylabel('r(footprint score, cold-ischaemia footprint score)\nper tumour cohort'); ax[1].set_ylim(-1, 1)
ax[1].set_title('B  Scores track the cold-ischaemia response', loc='left', fontsize=10)
for a in ax:
    a.set_facecolor('#fcfcfb'); [a.spines[s].set_visible(False) for s in ['top', 'right']]; a.grid(axis='y', color='#e5e4df', lw=.6); a.set_axisbelow(True)
plt.tight_layout(); plt.savefig(RA + 'fig_audit.png', dpi=150, facecolor='#fcfcfb')
