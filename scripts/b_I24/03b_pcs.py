"""Global phospho PCs per cohort (non-cascade sites, >=80% coverage, NaN->0 after z) for adjustment."""
from common import *
A = load_stages('tc_padj'); stage_prot = set(U.values()); out = []
for c in COHORTS:
    D = A[A.cohort == c]; X = cl.matrix(c, 'input_tc_padj')[D.index]
    X = X[[k.split('_')[0] not in stage_prot for k in X.index]]; X = X[X.notna().mean(1) >= 0.8]
    Z = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0).fillna(0).T.values
    U_, s, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    pc = pd.DataFrame(U_[:, :10] * s[:10], index=D.index, columns=[f'PC{i+1}' for i in range(10)])
    print(c, X.shape[0], 'var expl top5', (s[:5]**2 / (s**2).sum()).round(3))
    out.append(pc)
pd.concat(out).to_parquet(f'{OUT}/tables/phospho_pcs.parquet')
