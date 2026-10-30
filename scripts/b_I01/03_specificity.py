"""Family-echo quantification, specific residual, split-half reliability vs matched random sets, inhibitor positive control."""
import numpy as np, pandas as pd, scipy.sparse as sp, sys, os
from pathlib import Path
from scipy.stats import mannwhitneyu
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
TAG = os.environ.get("TAG", "")
SET_TAG = os.environ.get("SET_TAG", TAG)
PCS = int(os.environ.get("PCS", 0))  # if >0 also partial out top PCs of the kinome-wide full-score matrix
rng = np.random.default_rng(0)
NSPLIT, NNULL, NMIN = 10, 50, 5

conds = pd.read_csv(CA / "z_conds.tsv", sep="\t").iloc[:, 0].tolist()
C = load.conditions().set_index("conditionID").loc[conds]
study = C.study_pmid.values
ust, sidx = np.unique(study, return_inverse=True)
Zf = np.load(CA / "Zf.npy"); obs = np.load(CA / "obs.npy").astype(np.float32)
sets = np.load(CA / f"sets{SET_TAG}.npy", allow_pickle=True).item()
M = pd.read_csv(OUT / f"tables/substrate_sets{SET_TAG}.tsv", sep="\t", index_col=0)
K = load.kinases().set_index("kinase")
kin = list(M.index)


def score(idx_lists):
    rows = np.concatenate([np.asarray(ix, int) for ix in idx_lists])
    cols = np.concatenate([np.full(len(ix), j) for j, ix in enumerate(idx_lists)])
    I = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(Zf.shape[0], len(idx_lists)))
    s = np.asarray(I.T @ Zf); n = np.asarray(I.T @ obs)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = s / np.sqrt(n)
    out[n < NMIN] = np.nan
    return out


def scenter(A):
    """study-centre rows of A (sets x conds); drop studies with <3 non-NaN conds."""
    A = A.copy()
    for s in range(len(ust)):
        m = sidx == s
        sub = A[:, m]
        cnt = (~np.isnan(sub)).sum(1, keepdims=True)
        mu = np.nanmean(np.where(cnt >= 3, sub, np.nan), 1, keepdims=True) if m.sum() else 0
        sub = np.where(cnt >= 3, sub - mu, np.nan)
        A[:, m] = sub
    return A


XPC = None


def resid(y, x):
    """OLS residual of y on x (centred; plus global PCs if PCS>0), NaN-aware. Returns (residual, beta of x)."""
    m = ~np.isnan(y) & ~np.isnan(x)
    if m.sum() < 20:
        return y.copy(), np.nan
    X = x[:, None] if XPC is None else np.c_[x, XPC]
    b = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
    r = y - np.nan_to_num(X, nan=0.0) @ b
    return r, b[0]


def pcor(a, b):
    m = ~np.isnan(a) & ~np.isnan(b)
    return np.corrcoef(a[m], b[m])[0, 1] if m.sum() >= 20 else np.nan


import warnings; warnings.filterwarnings("ignore")
full = scenter(pd.read_parquet(CA / f"score_full{SET_TAG}.parquet").loc[kin].values)
spec = scenter(pd.read_parquet(CA / f"score_spec{SET_TAG}.parquet").loc[kin].values)
conf = scenter(pd.read_parquet(CA / f"score_conf{SET_TAG}.parquet").loc[kin].values)
# compendium core-v3 z_core (study-centred) for reference
L = load.activity_long(); L = L[(L.method == "corev3")]
zc = L.pivot_table(index="kinase", columns="condition", values="z_core").reindex(index=kin, columns=conds)
ns = L.pivot_table(index="kinase", columns="condition", values="n_substrates").reindex(index=kin, columns=conds)
zc_all = scenter(zc.values.astype(float))
zc = scenter(zc.where(ns >= 1).values.astype(float))
if PCS > 0:
    F0 = np.nan_to_num(full, nan=0.0); F0 = F0 - F0.mean(1, keepdims=True)
    U, Sv, Vt = np.linalg.svd(F0, full_matrices=False)
    XPC = Vt[:PCS].T  # conds x PCS
    print("PC var explained", (Sv[:PCS] ** 2 / (Sv ** 2).sum()).round(3), "total", (Sv[:PCS] ** 2).sum() / (Sv ** 2).sum())

# group-matched null pool: sites in top1% of any kinase of the same Manning group & residue class, excluding k's full set
P = pd.read_parquet(CA / "pct.parquet")
sites_z = pd.read_csv(CA / "z_sites.tsv", sep="\t").iloc[:, 0].values
P = P.loc[sites_z]
grp = K.loc[kin, "manning_group"].values
top = np.nan_to_num(P.values, nan=0) >= 99
pool_g = {g: np.where(top[:, grp == g].any(1))[0] for g in np.unique(grp)}
nonnan = ~np.isnan(P.values)
prot = np.array([x.split("_")[0] for x in sites_z])

rows = []; RES = np.full_like(spec, np.nan)
batch_sets, batch_meta = [], []
for i, k in enumerate(kin):
    sp_ix = sets[k]["spec"]
    r, b = resid(spec[i], conf[i]); RES[i] = r
    d = dict(kinase=k, n_spec=len(sp_ix), n_cond_spec=int((~np.isnan(spec[i])).sum()),
             r2_full_conf=pcor(full[i], conf[i]) ** 2, r2_spec_conf=pcor(spec[i], conf[i]) ** 2,
             r_spec_full=pcor(spec[i], full[i]), r2_corev3_conf=pcor(zc_all[i], conf[i]) ** 2,
             r_resid_corev3=pcor(r, zc_all[i]), beta_conf=b)
    rows.append(d)
    if len(sp_ix) < 20:
        continue
    excl = set(sets[k]["full"])
    pool = np.array([s for s in pool_g[grp[i]] if s not in excl and nonnan[s, i]])
    prots = prot[sp_ix]; up = np.unique(prots)
    for t in range(NSPLIT):  # split by PROTEIN so sites of one protein never sit in both halves
        a = set(rng.permutation(up)[: len(up) // 2]); ina = np.array([p in a for p in prots])
        batch_sets += [sp_ix[ina], sp_ix[~ina]]; batch_meta += [(i, "obs", t)] * 2
    for t in range(NNULL):
        if len(pool) < len(sp_ix): break
        rs = rng.choice(pool, len(sp_ix), replace=False); pr = prot[rs]; upr = np.unique(pr)
        a = set(rng.permutation(upr)[: len(upr) // 2]); ina = np.array([p in a for p in pr])
        batch_sets += [rs[ina], rs[~ina]]; batch_meta += [(i, "null", t)] * 2
print("scoring", len(batch_sets), "half-sets")
H = scenter(score(batch_sets))
bm = pd.DataFrame(batch_meta, columns=["i", "kind", "t"])
rel = {}
for (i, kind, t), g in bm.groupby(["i", "kind", "t"]):
    a, b2 = H[g.index[0]], H[g.index[1]]
    a, _ = resid(a, conf[i]); b2, _ = resid(b2, conf[i])
    rel.setdefault((i, kind), []).append(pcor(a, b2))
T = pd.DataFrame(rows).set_index("kinase")
for i, k in enumerate(kin):
    o = rel.get((i, "obs")); n = rel.get((i, "null"))
    if o is None: continue
    o = np.nanmean(o); n = np.array(n, float) if n else np.array([np.nan])
    T.loc[k, "rel_obs"] = o; T.loc[k, "rel_null_mean"] = np.nanmean(n); T.loc[k, "rel_null_sd"] = np.nanstd(n)
    T.loc[k, "rel_p"] = (1 + np.sum(n >= o)) / (1 + np.isfinite(n).sum())
    # Spearman-Brown full-length reliability
    T.loc[k, "rel_sb"] = 2 * o / (1 + o) if o > -1 else np.nan
T["rel_z"] = (T.rel_obs - T.rel_null_mean) / T.rel_null_sd
T["dark"] = K.loc[T.index, "idg_dark_kinase"].values
T["family"] = K.loc[T.index, "manning_family"].values
T["group"] = K.loc[T.index, "manning_group"].values
T.to_csv(OUT / f"tables/specificity_per_kinase{TAG}.tsv", sep="\t")
pd.DataFrame(RES, index=kin, columns=conds).to_parquet(CA / f"resid{TAG}.parquet")
pd.DataFrame(spec, index=kin, columns=conds).to_parquet(CA / f"spec_sc{TAG}.parquet")
pd.DataFrame(conf, index=kin, columns=conds).to_parquet(CA / f"conf_sc{TAG}.parquet")
pd.DataFrame(full, index=kin, columns=conds).to_parquet(CA / f"full_sc{TAG}.parquet")
pd.DataFrame(zc_all, index=kin, columns=conds).to_parquet(CA / f"corev3_sc.parquet")
print(T.groupby("dark")[["n_spec", "r2_full_conf", "r2_spec_conf", "r2_corev3_conf", "r_resid_corev3", "rel_obs", "rel_null_mean", "rel_z"]].median())
print("frac rel_p<0.05:", T.groupby("dark").apply(lambda d: ((d.rel_p <= 0.05)).mean()))

# ---- positive control: inhibitor primary / profiled targets
ok = C.target_confidence.isin(["high", "medium"]).values & (C.dose_nM.fillna(0).values >= 100) & (C.perturbation_class == "kinase_inhibitor").values & ~C.contrast_reversed.astype(bool).values
prim = C.inhibitor_primary_targets_in_universe.fillna("").str.split(";")
prof = C.inhibitor_profiled_targets.fillna("").str.split(";")
pc = []
for i, k in enumerate(kin):
    for lab, tg in [("primary", prim), ("profiled", prof)]:
        hit = ok & tg.apply(lambda x: k in x).values
        # exclude conditions where any confounder kinase is a primary target (to isolate k-specific effect)
        conf_k = set(M.loc[k, "neighbours"].split(";")) | set(K.index[(K.manning_family == K.loc[k, "manning_family"]) & (K.index != k)])
        conf_hit = prim.apply(lambda x: bool(conf_k & set(x))).values | prof.apply(lambda x: bool(conf_k & set(x))).values
        for clean in [False, True]:
            h = hit & (~conf_hit if clean else True)
            nst = len(set(study[h]))
            if h.sum() < 3 or nst < 1: continue
            d = dict(kinase=k, target_def=lab, conf_hit_excluded=clean, n_cond=int(h.sum()), n_study=nst)
            for nm, A in [("full", full), ("spec", spec), ("resid", RES), ("conf", conf), ("corev3", zc_all)]:
                a = A[i]; m = ~np.isnan(a)
                x1, x0 = a[h & m], a[(~hit) & m & ok]
                if len(x1) >= 3 and len(x0) >= 10:
                    d[f"auc_{nm}"] = 1 - mannwhitneyu(x1, x0).statistic / (len(x1) * len(x0))  # P(target lower)
                    d[f"n_{nm}"] = len(x1)
            pc.append(d)
pc = pd.DataFrame(pc); pc["dark"] = K.loc[pc.kinase, "idg_dark_kinase"].values
pc.to_csv(OUT / f"tables/inhibitor_control{TAG}.tsv", sep="\t", index=False)
print(pc.groupby(["target_def", "conf_hit_excluded", "dark"])[[c for c in pc.columns if c.startswith("auc_")]].median())
print(pc.groupby(["target_def", "conf_hit_excluded", "dark"]).size())
