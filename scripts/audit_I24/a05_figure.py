"""Audit figure: RTK-driven (IV) vs RTK-independent adaptor->ERK transmission across specifications."""
from common_a import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M1 = pd.read_csv(f'{OUT}/tables/A1_decomp_meta.tsv', sep='\t'); M2 = pd.read_csv(f'{OUT}/tables/A2_protein_dissection_meta.tsv', sep='\t')
M3 = pd.read_csv(f'{OUT}/tables/A3_abundance_split_meta.tsv', sep='\t'); G4 = pd.read_csv(f'{OUT}/tables/A4b_compendium_EGFR_transitivity.tsv', sep='\t')
def g(M, a, v='tc'): return M[(M.analysis == a) & (M.variant == v)].iloc[0]
rows = [('Unadjusted pY (I24 headline)', g(M1, 'main')),
        ('Unadjusted, + adaptor & ERK protein', g(M2, 'tc_adj_ADAPT_ERKprot')),
        ('Unadjusted, + RTK protein', g(M2, 'tc_adj_RTKprot')),
        ('Unadjusted, instrument = RTK pY/protein', g(M3, 'instr_RTK_stoich')),
        ('Protein-adjusted pY', g(M1, 'main', 'tc_padj')),
        ('Protein-adjusted, EGFR pY only', g(M1, 'EGFR_only', 'tc_padj'))]
c1, c2, ink, muted = '#2a78d6', '#eb6834', '#1f1f1e', '#6b6a63'
fig, ax = plt.subplots(figsize=(9, 4.8), dpi=150)
for i, (lab, r) in enumerate(rows):
    y = len(rows) - i
    ax.errorbar(r.iv, y + 0.13, xerr=[[r.iv - r.iv_lo], [r.iv_hi - r.iv]], fmt='o', color=c1, ms=6, lw=2, capsize=0)
    ax.errorbar(r.other, y - 0.13, xerr=[[r.other - r.other_lo], [r.other_hi - r.other]], fmt='s', color=c2, ms=6, lw=2, capsize=0)
    ax.text(1.02, y, f'p={r.diff_p:.2g}', transform=ax.get_yaxis_transform(), va='center', fontsize=8, color=muted)
e = G4.set_index('lab').loc['EGFRaxis_all']
ax.errorbar(e.iv, 0.13, xerr=[[e.iv - e.iv_lo], [e.iv_hi - e.iv]], fmt='o', color=c1, ms=6, lw=2)
ax.plot(e.other, -0.13, 's', color=c2, ms=6)
ax.text(1.02, 0, 'n=47 cond.', transform=ax.get_yaxis_transform(), va='center', fontsize=8, color=muted)
ax.set_yticks(range(0, len(rows) + 1)); ax.set_yticklabels(['Cell lines: EGF / EGFR-i (compendium)'] + [r[0] for r in rows][::-1], fontsize=8.5, color=ink)
ax.axvline(0, color='#c9c8c0', lw=1); ax.set_xlim(-0.6, 1.3)
ax.set_xlabel('Adaptor pY → ERK transmission (SD/SD), meta over 7 CPTAC cohorts, 95% CI', fontsize=8.5, color=ink)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.plot([], [], 'o', color=c1, label='RTK-driven part (2SLS, instrument RTK pY)'); ax.plot([], [], 's', color=c2, label='RTK-independent part (residual)')
ax.legend(fontsize=8, frameon=False, loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=2)
ax.set_title('Non-transitivity gap is confined to unadjusted pY', fontsize=10, color=ink, loc='left', pad=22)
plt.tight_layout(rect=(0,0,0.95,1)); plt.savefig(f'{OUT}/figures/AF1_decomposition_specifications.png'); print('ok')
