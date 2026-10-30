"""Kinase-kinase perturbation association matrices (query q x target j), across compendium conditions.
a_resid : corr(resid_q, score of full_j with q's specific sites removed)      [kinase-specific, family partialled]
a_full  : corr(full_q, full_j)                                                [naive, standard activity]
a_conf  : corr(conf_q, full_j)                                                [family/neighbour-only baseline]
a_corev3: corr(core-v3 z_q, core-v3 z_j)                                      [compendium standard]
Correlation: study-centred values, NaN->0, condition weight 1/sqrt(n conditions in study) (limits decryptM dominance)."""
import numpy as np, pandas as pd, scipy.sparse as sp, sys, os
from pathlib import Path
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
TAG = os.environ.get("TAG", "")
SET_TAG = os.environ.get("SET_TAG", TAG)
conds = pd.read_csv(CA / "z_conds.tsv", sep="\t").iloc[:, 0].tolist()
C = load.conditions().set_index("conditionID").loc[conds]
study = C.study_pmid.values
ust, sidx, cnt = np.unique(study, return_inverse=True, return_counts=True)
w = 1 / np.sqrt(cnt[sidx]); w = w / w.sum()
Zf = np.load(CA / "Zf.npy"); obs = np.load(CA / "obs.npy").astype(np.float32)
sets = np.load(CA / f"sets{SET_TAG}.npy", allow_pickle=True).item()
kin = list(pd.read_csv(OUT / f"tables/substrate_sets{SET_TAG}.tsv", sep="\t", index_col=0).index)
nk = len(kin)


def scenter(A):
    A = A.copy()
    for s in range(len(ust)):
        m = sidx == s; sub = A[:, m]
        c = (~np.isnan(sub)).sum(1, keepdims=True)
        with np.errstate(all="ignore"):
            mu = np.nanmean(np.where(c >= 3, sub, np.nan), 1, keepdims=True)
        A[:, m] = np.where(c >= 3, sub - mu, np.nan)
    return A


def wstd(A):
    """weighted standardisation of rows (NaN->0 after centring)"""
    A = np.nan_to_num(A, nan=0.0)
    mu = (A * w).sum(1, keepdims=True); A = A - mu
    sd = np.sqrt((A ** 2 * w).sum(1, keepdims=True)); sd[sd == 0] = np.nan
    return A / sd


def wcorr(A, B):
    return (wstd(A) * w) @ wstd(B).T


# full_j raw sums & counts
def IM(lists):
    rows = np.concatenate([np.asarray(x, int) for x in lists]) if lists else np.array([], int)
    cols = np.concatenate([np.full(len(x), j) for j, x in enumerate(lists)])
    return sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(Zf.shape[0], len(lists)))


If = IM([sets[k]["full"] for k in kin])
Sf = np.asarray(If.T @ Zf); Nf = np.asarray(If.T @ obs)
res = pd.read_parquet(CA / f"resid{TAG}.parquet").loc[kin].values
full = pd.read_parquet(CA / f"full_sc{SET_TAG}.parquet").loc[kin].values
conf = pd.read_parquet(CA / f"conf_sc{SET_TAG}.parquet").loc[kin].values
cv3 = pd.read_parquet(CA / "corev3_sc.parquet").loc[kin].values
with np.errstate(all="ignore"):
    fullJ = np.where(Nf >= 5, Sf / np.sqrt(Nf), np.nan)
fullJs = wstd(scenter(fullJ))
A_full = wcorr(full, full); A_conf = wcorr(conf, full); A_cv3 = wcorr(cv3, cv3)
A_res = np.full((nk, nk), np.nan)
Rs = wstd(res)
for qi, q in enumerate(kin):
    spq = set(sets[q]["spec"])
    if len(spq) < 20 or np.all(np.isnan(res[qi])): continue
    inter = [np.array([s for s in sets[j]["full"] if s in spq], int) for j in kin]
    Ii = IM(inter)
    Si = np.asarray(Ii.T @ Zf); Ni = np.asarray(Ii.T @ obs)
    S2, N2 = Sf - Si, Nf - Ni
    with np.errstate(all="ignore"):
        J = np.where(N2 >= 5, S2 / np.sqrt(N2), np.nan)
    Js = wstd(scenter(J))
    A_res[qi] = (Rs[qi] * w) @ Js.T
    if qi % 50 == 0: print(qi, flush=True)
for nm, A in [("resid", A_res), ("full", A_full), ("conf", A_conf), ("corev3", A_cv3)]:
    pd.DataFrame(A, index=kin, columns=kin).to_parquet(CA / f"assoc_{nm}{TAG}.parquet")
