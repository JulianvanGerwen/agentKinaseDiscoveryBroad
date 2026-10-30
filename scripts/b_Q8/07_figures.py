"""b_Q8 figures F1-F4."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b'})
BLUE, ORANGE, AQUA, GREY = '#2a78d6', '#eb6834', '#1baf7a', '#9a9893'
CL = {'GAINED in tumour': ORANGE, 'LOST/attenuated in tumour': BLUE, 'preserved': AQUA}
P = pd.read_csv(TQ + 'T2_coherence_classes_primary.tsv', sep='\t', index_col=0)
C = pd.read_csv(TQ + 'T1_coherence.tsv', sep='\t')
# F1
fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
P['cls2'] = P['class'].where(P['class'].isin(list(CL)), 'no robust Δ / incoherent')
for cl_ in ['no robust Δ / incoherent', 'preserved', 'GAINED in tumour', 'LOST/attenuated in tumour']:
    d = P[P.cls2 == cl_]
    if len(d): ax[0].scatter(d.z_N_med, d.z_T_med, s=22, c=CL.get(cl_, GREY), label=cl_, edgecolor='white', lw=0.6)
lim = [min(P.z_N_med.min(), P.z_T_med.min()) - 1, max(P.z_N_med.max(), P.z_T_med.max()) + 1]
ax[0].plot(lim, lim, c=GREY, lw=0.8, ls='--'); ax[0].set_xlabel('NAT coherence z vs matched random sites (median, 6 cohorts)')
ax[0].set_ylabel('Tumour coherence z (paired tumours, same n)')
for s, r in P.iterrows():
    if r['class'] == 'GAINED in tumour' or max(r.z_N_med, r.z_T_med) > 8: ax[0].annotate(s.split(':')[1], (r.z_N_med, r.z_T_med), fontsize=6.5, color='#52514e', xytext=(3, 2), textcoords='offset points')
ax[0].legend(frameon=False, fontsize=7); ax[0].set_title('Footprint coherence, NAT vs tumour (protein- + 3-PC-adjusted)', loc='left')
d = C[(C.cond == 'padj_pc3') & C.cohort.isin(MAINC)].pivot(index='set', columns='cohort', values='z_delta_perm')
d = d.loc[P.index.intersection(d.index)]; d = d.loc[P.loc[d.index, 'meta_z_perm'].sort_values().index]
sel = d.index[(P.loc[d.index, 'meta_z_perm'].abs() > 2) | (P.loc[d.index, 'nN_sig'] >= 4)]
d = d.loc[sel]
im = ax[1].imshow(d.values, cmap='RdBu_r', vmin=-4, vmax=4, aspect='auto')
ax[1].set_yticks(range(len(d))); ax[1].set_yticklabels([s.split(':')[1] + (' (PSP)' if s.startswith('PSP') else '') for s in d.index], fontsize=6.5)
ax[1].set_xticks(range(d.shape[1])); ax[1].set_xticklabels(d.columns, rotation=45)
plt.colorbar(im, ax=ax[1], label='Δ coherence (T−NAT), permutation z', shrink=0.7)
ax[1].set_title('Per-cohort differential coherence', loc='left')
plt.tight_layout(); plt.savefig(FQ + 'F1_coherence_NAT_vs_tumour.png', dpi=200); plt.close()
# F2 edges
E = pd.read_csv(TQ + 'T3_compendium_edges.tsv', sep='\t'); E['rel'] = E.a + '|' + E.b; E = E.set_index('rel')
D = pd.read_parquet(CQ + 'edge_stats_padj.parquet'); D = D[D.kind == 'pair']; D['label'] = D.rel.map(E.label); D['s'] = np.sign(D.rel.map(E.r_meta))
fig, axs = plt.subplots(2, 4, figsize=(11, 5.6), sharex=True, sharey=True)
for a, coh in zip(axs.ravel(), NATC):
    x = D[(D.cohort == coh)]
    ne = x[x.label == 'nonedge']; ed = x[x.label == 'edge']
    a.scatter(ne.r_N, ne.r_T, s=4, c=GREY, alpha=0.4, lw=0, label='compendium non-edge (|r|<.05)')
    a.scatter(ed.s * ed.r_N, ed.s * ed.r_T, s=6, c=BLUE, alpha=0.6, lw=0, label='compendium edge (sign-oriented)')
    a.axhline(0, c=GREY, lw=0.5); a.axvline(0, c=GREY, lw=0.5); a.plot([-1, 1], [-1, 1], c=GREY, lw=0.6, ls='--')
    a.set_title(f'{coh} (n={int(x.n.iloc[0])}/group)', loc='left', fontsize=8)
    a.text(0.03, 0.9, f'sign-conc N {np.mean(np.sign(ed.r_N) == ed.s):.2f} / T {np.mean(np.sign(ed.r_T) == ed.s):.2f}', transform=a.transAxes, fontsize=6.5)
for a in axs[1]: a.set_xlabel('r in NATs')
for a in axs[:, 0]: a.set_ylabel('r in tumours')
axs[0, 0].legend(frameon=False, fontsize=6, loc='lower right')
fig.suptitle('Compendium programme–programme relationships are preserved in both normal tissue and tumours (protein-adjusted)', x=0.01, ha='left')
plt.tight_layout(); plt.savefig(FQ + 'F2_edge_preservation.png', dpi=200); plt.close()
# F3 axis
A = pd.read_csv(TQ + 'T4_axis.tsv', sep='\t'); A = A[A.cond == 'padj']
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
for i, (col, ttl) in enumerate([('r_MEK_EGF_exclusive', 'MEKi vs EGF footprint r (no shared sites/proteins)'), ('PC1_var', 'PC1 variance of 49 footprint scores')]):
    p = A.pivot(index='cohort', columns='grp', values=col).loc[NATC]
    xx = np.arange(len(p)); ax[i].bar(xx - 0.2, p['N'], 0.38, color=AQUA, label='NAT'); ax[i].bar(xx + 0.2, p['T'], 0.38, color=ORANGE, label='tumour (paired, same n)')
    ax[i].set_xticks(xx); ax[i].set_xticklabels(p.index, rotation=45); ax[i].set_title(ttl, loc='left'); ax[i].axhline(0, c=GREY, lw=0.6)
ax[0].legend(frameon=False, fontsize=7)
plt.tight_layout(); plt.savefig(FQ + 'F3_MEK_EGF_axis.png', dpi=200); plt.close()
print('figures done')
