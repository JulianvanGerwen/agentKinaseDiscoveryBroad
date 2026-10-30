"""I27 step 5: figures."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T = '/workspace/results/b_I27/tables/'
F = '/workspace/results/b_I27/figures/'
R = pd.read_csv(T + 'T1_confirmatory.tsv', sep='\t')
FU = pd.read_csv(T + 'T6_followup.tsv', sep='\t')
E = pd.read_csv(T + 'T4_exploratory.tsv', sep='\t')

# F1: forest of pre-registered primary outcomes
panels = [
    ('P1 KRAS G12R vs G12D/V (PDAC)', R[(R.pred.isin(['P1', 'P1_primary', 'P1_matched'])) &
      R.y.isin(['stage_AKTSUB', 'stage_MTORC1', 'stage_S6', 'site_PRAS40_T246', 'kin_AKT1', 'stage_ERK', 'stage_ERKSUB'])]),
    ('P2 KRAS G12C vs G12D (LUAD)', R[(R.pred == 'P2_G12CvG12D') & (R.cohort == 'LUAD') &
      R.y.isin(['stage_AKTSUB', 'stage_MTORC1', 'site_PRAS40_T246', 'kin_AKT1', 'stage_ERK', 'stage_ERKSUB'])]),
    ('P3 PIK3CA helical vs kinase (meta)', None),
    ('P4 BRAF V600 vs RAS-mut (COAD)', R[(R.pred == 'P4_V600_v_RASmut') &
      R.y.isin(['stage_RTK', 'stage_ADAPT', 'stage_EGFRpY', 'stage_ERK', 'stage_ERKSUB', 'stage_AKTSUB'])]),
    ('P5 EGFR ex19del vs L858R (LUAD)', R[(R.pred == 'P5_ex19del_v_L858R') &
      R.y.isin(['stage_RTK', 'stage_EGFRpY', 'stage_ADAPT', 'stage_ERK', 'stage_ERKSUB', 'stage_AKTSUB', 'stage_MTORC1'])]),
]
try:
    M = pd.read_csv(T + 'T3_meta.tsv', sep='\t')
    mm = M[(M.pred == 'P3_helical_v_kinase') &
           M.y.isin(['SGKmAKT', 'stage_AKTSUB', 'stage_MTORC1', 'site_PRAS40_T246', 'site_NDRG1_SGK', 'kin_SGK3'])].copy()
    mm = mm.rename(columns={'beta_re': 'beta', 'lo_re': 'lo', 'hi_re': 'hi', 'p_re': 'p'})
    panels[2] = (panels[2][0], mm)
except Exception:
    panels[2] = (panels[2][0], pd.DataFrame())

fig, axes = plt.subplots(1, 5, figsize=(19, 5.2))
for ax, (title, d) in zip(axes, panels):
    if d is None or not len(d):
        ax.set_title(title + '\n(no estimable)'); ax.axis('off'); continue
    d = d.dropna(subset=['beta']).drop_duplicates('y').sort_values('y')
    yv = np.arange(len(d))
    ax.errorbar(d.beta, yv, xerr=[d.beta - d.lo, d.hi - d.beta], fmt='o', color='#2b6cb0', ms=5, lw=1.4, capsize=3)
    ax.axvline(0, color='#888', lw=1, ls='--')
    ax.set_yticks(yv); ax.set_yticklabels([y.replace('stage_', '').replace('site_', '').replace('kin_', '') for y in d.y], fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('effect (within-cohort SD)')
    for k, (_, r) in enumerate(d.iterrows()):
        if pd.notna(r.p) and r.p < 0.05:
            ax.text(r.hi + 0.05, k, '*', va='center', fontsize=12, color='#c53030')
    ax.set_xlim(-3, 3)
fig.suptitle('I27: pre-registered allele contrasts (HC3 95% CI; * permutation p<.05)', fontsize=12)
fig.tight_layout()
fig.savefig(F + 'F1_prereg_forest.png', dpi=140)

# F2: LUAD MAPK activation-loop axis by allele, with smoking adjustment
d = FU[(FU.cohort == 'LUAD') & FU.y.isin(['kin_MAP2Kfam', 'kin_MEK12', 'loop_ERK_TEY', 'stage_ERKSUB', 'fp_KI_MEK'])]
fig, ax = plt.subplots(figsize=(9, 5))
labs = []
k = 0
for (a, b), g in d.groupby(['a', 'b'], sort=False):
    for y in ['kin_MAP2Kfam', 'kin_MEK12', 'loop_ERK_TEY', 'stage_ERKSUB', 'fp_KI_MEK']:
        for mdl, off, col in [('base', -0.15, '#2b6cb0'), ('smoke_adj', 0.15, '#dd6b20')]:
            r = g[(g.y == y) & (g.model == mdl)]
            if not len(r): continue
            r = r.iloc[0]
            ax.errorbar(r.beta, k + off, xerr=[[r.beta - r.lo], [r.hi - r.beta]], fmt='o', color=col, ms=4, lw=1.2, capsize=2)
        labs.append(f'{a} vs {b}: {y}')
        k += 1
ax.axvline(0, color='#888', ls='--', lw=1)
ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=7.5)
ax.set_xlabel('effect (within-cohort SD)')
ax.set_title('LUAD: MAPK-module outputs by driver allele (blue = base, orange = smoking-adjusted)', fontsize=10)
fig.tight_layout(); fig.savefig(F + 'F2_luad_mapk.png', dpi=140)

# F3: exploratory volcano
fig, ax = plt.subplots(figsize=(7.5, 5.5))
col = np.where(E.q < 0.1, '#c53030', '#aaa')
ax.scatter(E.beta, -np.log10(E.p), s=8, c=col, alpha=.6, lw=0)
ax.set_xlabel('effect (within-cohort SD)'); ax.set_ylabel('-log10 p')
ax.set_title(f'Exploratory allele contrasts (n={len(E)} tests, {int((E.q<.1).sum())} at FDR 10%)', fontsize=10)
for _, r in E.nsmallest(8, 'p').iterrows():
    ax.annotate(f"{r.cohort} {r.a}/{r.b}\n{r.y}", (r.beta, -np.log10(r.p)), fontsize=6.5,
                xytext=(4, 2), textcoords='offset points')
fig.tight_layout(); fig.savefig(F + 'F3_exploratory_volcano.png', dpi=140)
print('figures done')
