"""audit_I10 step 6: is 'PC1 = 27%' of the 49 footprint scores a property of the footprints, or of the scoring/sign pattern?
Uses b_I10's own cached scores (tc_padj nz) and its cached matched-random-set raw scores (cptac_null_*.npy, same weights/signs, 200 draws).
Null PC1: for draw k, every footprint replaced by its k-th random matched set; nz computed exactly as b_I10 (vs the other 199 draws);
within-cohort z; PCA; PC1 variance. Also: PC1 loading vs fraction of down-weighted sites; PC1 vs MAPK1/3; PC1 after removing MAPK-redundant sets.
Output: A6_pc1_null.tsv"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from common import *
names = pd.read_csv('/workspace/cache/b_I10/cptac_score_names.tsv', sep='\t', header=None)[0].tolist()
S = pd.read_parquet('/workspace/cache/b_I10/cptac_scores_tc_padj.parquet'); S = S[~S.footprint.str.contains(r'\|')]
fps = sorted(S.footprint.unique()); jj = [names.index(f) for f in fps]
fd = FP.groupby('group').w.apply(lambda w: (w < 0).mean())
def pca(Fz):
    U, s, Vt = np.linalg.svd(Fz, full_matrices=False); return (s ** 2 / (s ** 2).sum()), Vt[0], U[:, 0] * s[0]
blocks, nblocks = [], {k: [] for k in range(0, 200, 10)}
for coh in COH:
    tum = pd.read_parquet(f'/workspace/cache/b_I10/cptac_rank_tc_padj_{coh}.parquet').columns
    F = S[S.cohort == coh].pivot(index='sample', columns='footprint', values='nz').reindex(tum)[fps]
    blocks.append(((F - F.mean()) / F.std()).fillna(0))
    N = np.load(f'/workspace/cache/b_I10/cptac_null_tc_padj_{coh}.npy')[jj]          # fp x NR x samples
    for k in nblocks:
        others = np.delete(N, k, axis=1)
        nz = (N[:, k] - np.nanmean(others, 1)) / np.nanstd(others, 1)              # fp x samples
        D = pd.DataFrame(nz.T, index=tum, columns=fps); nblocks[k].append(((D - D.mean()) / D.std()).fillna(0))
Fz = pd.concat(blocks); ve, load_, pc1 = pca(Fz.values)
L = pd.Series(load_, index=fps)
if L['KI_MEK'] > 0: L, pc1 = -L, -pc1
nve = np.array([pca(pd.concat(b).values)[0][0] for b in nblocks.values()])
mk = cl.activity('tc_padj', 'corev3_z').loc[['MAPK1', 'MAPK3']].mean(0).reindex(Fz.index)
rows = [dict(what='PC1 variance, real footprints (b_I10 nz)', value=ve[0]),
        dict(what='PC1 variance, matched random sets with identical weights/signs (20 draws) median [min-max]', value=np.median(nve), note=f'{nve.min():.3f}-{nve.max():.3f}'),
        dict(what='Spearman(PC1 loading, fraction of down-weighted sites) across 49 footprints', value=spearmanr(L, fd.reindex(fps))[0]),
        dict(what='Spearman(PC1 sample score, MAPK1/3 corev3 tc_padj)', value=spearmanr(pc1, mk, nan_policy='omit')[0])]
mapk = ['EGF', 'GF_all', 'KI_EGFR', 'KI_MET_ALK', 'KI_RAF', 'KI_p38', 'MAPKaxis', 'KI_ERK', 'FGF', 'serum']
keep = [f for f in fps if f not in mapk]
v2, l2, p2 = pca(Fz[keep].values); L2 = pd.Series(l2, index=keep)
rows.append(dict(what='PC1 variance after dropping MAPK-redundant footprints (keep KI_MEK)', value=v2[0], note=';'.join(f'{k}:{v:+.2f}' for k, v in L2.sort_values(key=abs, ascending=False)[:10].items())))
keep2 = [f for f in keep if f != 'KI_MEK']
v3, l3, p3 = pca(Fz[keep2].values)
rows.append(dict(what='PC1 variance with no MAPK footprint at all', value=v3[0], note=f'rho(PC1, MAPK1/3)={abs(spearmanr(p3, mk, nan_policy="omit")[0]):.2f}; rho(loading, frac_down)={spearmanr(l3, fd.reindex(keep2))[0]:.2f}'))
T = pd.DataFrame(rows); T.to_csv(RA + 'A6_pc1_null.tsv', sep='\t', index=False)
pd.set_option('display.max_colwidth', 250); pd.set_option('display.width', 250); print(T.to_string())
print(pd.DataFrame({'load': L, 'frac_down': fd.reindex(fps)}).sort_values('load').round(2).T.to_string())
