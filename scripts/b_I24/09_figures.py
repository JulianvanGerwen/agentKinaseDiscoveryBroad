from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B, O, A_, G = '#2a78d6', '#eb6834', '#1baf7a', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': G, 'axes.labelcolor': '#0b0b0b'})
T = f'{OUT}/tables'
M = pd.read_csv(f'{T}/T2_gain_meta.tsv', sep='\t'); M = M[M.variant == 'tc_padj']
IV = pd.read_csv(f'{T}/T6_reliability_iv.tsv', sep='\t').groupby('edge')[['r_disatt', 'iv_halfA_by_B', 'iv_halfB_by_A']].median()
NL = pd.read_csv(f'{T}/T3_random_site_null_pc0.tsv', sep='\t').groupby('edge').null_mean.median()
order = ['RTK>ADAPT', 'ADAPT>ERK', 'ERK>RSK', 'ERK>ERKSUB', 'RTK>ERK', 'ADAPT>AKTSUB', 'AKTSUB>MTORC1', 'MTORC1>S6', 'JNK>JNKSUB', 'P38>P38SUB']
M = M.set_index('edge').loc[order]
fig, ax = plt.subplots(figsize=(6.5, 4)); y = np.arange(len(order))[::-1]
ax.errorbar(M.std_slope_meta, y, xerr=[M.std_slope_meta - M.std_lo, M.std_hi - M.std_slope_meta], fmt='o', color=B, ms=6, lw=2, label='observed gain (RE meta, 95% CI)')
iv = IV.reindex(order)[['iv_halfA_by_B', 'iv_halfB_by_A']].mean(1)
ax.plot(iv, y, 'D', color=O, ms=6, label='IV (split-half instrument), median')
ax.plot(NL.reindex(order), y, '|', color=G, ms=12, mew=2, label='random-site null, median')
ax.set_yticks(y); ax.set_yticklabels(order); ax.axvline(0, color=G, lw=0.8)
ax.set_xlabel('standardized gain (SD output per SD input; adj. purity, proliferation)')
ax.set_title('Stage-to-stage transmission gain, 10 CPTAC cohorts (protein-adjusted sites)', fontsize=9, loc='left')
ax.legend(frameon=False, fontsize=8, loc='lower right'); fig.tight_layout(); fig.savefig(f'{OUT}/figures/F1_gain_by_edge.png', dpi=160)

# F2: primary modifiers - interaction vs level
D = pd.read_csv(f'{T}/T4_modifiers_meta.tsv', sep='\t'); D = D[(D.model == 'main') & D.primary & (D.k >= 3)]
fig, ax = plt.subplots(1, 2, figsize=(9, 3.8), sharey=False)
for a, col, lo, hi, ttl in [(ax[0], 'int_meta', 'int_lo', 'int_hi', 'gain change (x:lesion)'), (ax[1], 'lev_meta', 'lev_lo', 'lev_hi', 'offset (output at fixed input)')]:
    d = D.sort_values(col).reset_index(drop=True); x = np.arange(len(d))
    sig = (d[col.replace('meta', 'q').replace('int_q', 'int_q_primary').replace('lev_q', 'lev_q_primary')] < 0.05) if True else None
    qcol = 'int_q_primary' if col.startswith('int') else 'lev_q_primary'; sig = d[qcol] < 0.05
    a.vlines(x, d[lo], d[hi], color=np.where(sig, O, '#b5b4ae'), lw=1)
    a.scatter(x, d[col], s=10, color=np.where(sig, O, B), zorder=3)
    a.axhline(0, color=G, lw=0.8); a.set_xticks([]); a.set_xlabel(f'{len(d)} pre-declared (edge, regulator, data type) tests, ranked')
    a.set_title(ttl + f' — {int(sig.sum())} at FDR<0.05', fontsize=9, loc='left'); a.set_ylim(-1.2, 1.2)
    if sig.any():
        txt = '\n'.join(f"{d.gene[i]} {d.type[i]} ({d.edge[i]}): {d[col][i]:+.2f}" for i in np.where(sig)[0])
        a.text(0.02, 0.97, 'FDR<0.05 (orange):\n' + txt, transform=a.transAxes, va='top', fontsize=7, color='#0b0b0b')
ax[0].set_ylabel('standardized effect (RE meta, 95% CI)')
fig.tight_layout(); fig.savefig(f'{OUT}/figures/F2_modifiers_gain_vs_offset.png', dpi=160)

# F3: GBM EGFR-amp ADAPT>ERK scatter + decomposition
A = load_stages('tc_padj'); Dg = A[A.cohort == 'GBM'].copy()
gi = cl.matrix('GBM', 'cnv_gistic'); mb = cl.matrix('GBM', 'mut_binary')
amp = gi.loc['EGFR'].reindex(Dg.case_id).values == 2; nf1 = mb.loc['NF1'].reindex(Dg.case_id).fillna(0).values == 1
fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
for msk, col, lab in [(~amp & ~nf1, G, 'other'), (amp, B, 'EGFR amp'), (nf1 & ~amp, O, 'NF1 mut')]:
    q = Dg[msk]; ax[0].scatter(q.ADAPT, q.ERK, s=14, color=col, label=f'{lab} (n={msk.sum()})', alpha=0.85, edgecolor='white', linewidth=0.5)
    if msk.sum() >= 10:
        b = np.polyfit(q.ADAPT, q.ERK, 1); xx = np.linspace(q.ADAPT.min(), q.ADAPT.max(), 10); ax[0].plot(xx, np.polyval(b, xx), color=col, lw=2)
ax[0].set_xlabel('adaptor pY (SHC1/GAB1/PTPN11), log2'); ax[0].set_ylabel('ERK1/2 pT/pY, log2'); ax[0].legend(frameon=False, fontsize=8)
ax[0].set_title('GBM: EGFR-amp tumours sit below the line (offset -0.44 SD, p=0.03)', fontsize=9, loc='left')
P = pd.read_csv(f'{T}/T7_adaptor_decomposition_meta.tsv', sep='\t')
lab = []; i = 0
for v in ['tc', 'tc_padj']:
    for out in ['ERK', 'AKTSUB', 'ERKSUB']:
        r = P[(P.variant == v) & (P.out == out)].iloc[0]
        ax[1].errorbar([i - 0.15], [r.rtk_part], yerr=[[r.rtk_part - r.rtk_lo], [r.rtk_hi - r.rtk_part]], fmt='o', color=B, lw=2, label='RTK-driven adaptor pY (2SLS)' if i == 0 else None)
        ax[1].errorbar([i + 0.15], [r.other_part], yerr=[[r.other_part - r.oth_lo], [r.oth_hi - r.other_part]], fmt='s', color=O, lw=2, label='RTK-independent adaptor pY' if i == 0 else None)
        lab.append(f'{out}\n({v})'); i += 1
ax[1].set_xticks(range(i)); ax[1].set_xticklabels(lab, fontsize=7.5); ax[1].axhline(0, color=G, lw=0.8)
ax[1].set_ylabel('output SD per adaptor SD (RE meta, 7 cohorts)'); ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title('RTK-driven adaptor pY is transmitted less', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/figures/F3_decoupling.png', dpi=160)
