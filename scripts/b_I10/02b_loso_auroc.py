"""I10 step 2b: AUROC summary of LOSO footprint scores (from cache/b_I10/loso_scores.parquet)."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
comp = {'MAPKaxis': ['KI_MEK', 'KI_ERK', 'KI_RAF'], 'PI3Kaxis': ['KI_PI3K', 'KI_PI3K_mTOR', 'KI_mTORcat', 'KI_AKT', 'KI_PDPK1', 'KI_rapalog', 'KI_S6K'],
        'GF_all': ['EGF', 'FGF', 'insulin', 'serum', 'starvation']}
R = pd.read_parquet(CA + 'loso_scores.parquet')
R.loc[R.exp_group == 'starvation', 'score'] *= np.where(R.loc[R.exp_group == 'starvation', 'footprint'] == 'GF_all', -1, 1)
ax = lambda g: next((c for c, gs in comp.items() if g in gs or g == c), g)
out = []
for (fpg, v), d in R.dropna(subset=['score']).groupby(['footprint', 'variant']):
    pos = comp.get(fpg, [fpg]); y = d.exp_group.isin(pos)
    if y.sum() < 2 or (~y).sum() < 5: continue
    d2 = d[y | (d.exp_group.map(ax) != ax(fpg))]; y2 = d2.exp_group.isin(pos)
    out.append(dict(footprint=fpg, variant=v, n_pos=int(y.sum()), n_neg=int((~y).sum()), n_heldout_studies=d[y].heldout_study.nunique(),
                    auroc=roc_auc_score(y, d.score), auroc_vs_other_axis=roc_auc_score(y2, d2.score)))
A = pd.DataFrame(out)
P = A.pivot(index='footprint', columns='variant', values='auroc_vs_other_axis').add_prefix('aurocOA_')
P = P.join(A.pivot(index='footprint', columns='variant', values='auroc').add_prefix('auroc_'))
P = P.join(A[A.variant == 'all'].set_index('footprint')[['n_pos', 'n_neg', 'n_heldout_studies']])
P.to_csv(RT + 'T2_loso_auroc.tsv', sep='\t'); print(P.round(2).sort_values('aurocOA_all', ascending=False).to_string())
print('median AUROC (vs other axis): all/down/up', P[['aurocOA_all', 'aurocOA_down', 'aurocOA_up']].median().round(3).to_dict())
ki = P.index.str.startswith('KI_') | P.index.str.endswith('axis')
print('inhibitor footprints: down vs up paired', P.loc[ki, ['aurocOA_down', 'aurocOA_up']].median().round(3).to_dict(), (P.loc[ki,'aurocOA_down'] > P.loc[ki,'aurocOA_up']).mean())
# rank of own footprint among footprints for each held-out exp
a = R[R.variant == 'all'].dropna(subset=['score'])
a = a[~a.footprint.isin(comp)]
a['z'] = a.groupby('footprint').score.transform(lambda s: (s - s.mean()) / s.std())
own = a[a.footprint == a.exp_group]
rk = a.groupby('exp').z.rank(ascending=False)
own_rank = rk[own.index]; nfp = a.groupby('exp').size()[own.exp].values
print('own-footprint rank among', np.median(nfp), 'footprints: median', own_rank.median(), 'top1', (own_rank == 1).mean().round(3), 'top3', (own_rank <= 3).mean().round(3), 'expected top3', (3 / nfp).mean().round(3))
