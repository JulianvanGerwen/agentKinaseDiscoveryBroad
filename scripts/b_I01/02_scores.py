"""Kinase-specific substrate sets and per-condition scores.
For each kinase k:
  full_k : sites with pct_k >= TOP (top (100-TOP)% of U_obs by kinaseSubMain release score)
  conf set C_k : Manning family members + 5 nearest 'motif neighbours' (Pearson of pct vectors over sites)
  spec_k : full_k sites where max pct over C_k < MARG  (k-preferred sites)
  conf_k : sites in top of any C_k member with pct_k < MARG (confounders' sites, disjoint from k)
Score = Stouffer of per-condition robust z (median/MAD), n >= NMIN measured sites.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, sys, json
from pathlib import Path
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
import os
TOP, MARG, NMIN, NNEI = float(os.environ.get("TOP", 99)), float(os.environ.get("MARG", 99)), 5, 5
TAG = os.environ.get("TAG", "")

C = load.conditions(); C = C[C.qc_flag == "ok"]
K = load.kinases().set_index("kinase")
P = pd.read_parquet(CA / "pct.parquet")
X = load.sites_matrix(C.conditionID.tolist(), sites=P.index)
P = P.loc[X.index]
print("sites", X.shape)
# robust z per condition
Xv = X.values.astype(np.float32)
med = np.nanmedian(Xv, 0); mad = np.nanmedian(np.abs(Xv - med), 0) * 1.4826
Z = np.clip((Xv - med) / mad, -6, 6)
obs = ~np.isnan(Z); Zf = np.where(obs, Z, 0).astype(np.float32)
np.save(CA / "Zf.npy", Zf); np.save(CA / "obs.npy", obs)
pd.Series(X.index).to_csv(CA / "z_sites.tsv", index=False, sep="\t"); pd.Series(X.columns).to_csv(CA / "z_conds.tsv", index=False, sep="\t")

kin = list(P.columns)
Pv = P.values
# motif neighbours: correlation of pct vectors (sites with both non-NaN; use ST/Y separately implicitly via NaN->mean)
Pm = np.where(np.isnan(Pv), np.nanmean(Pv, 0), Pv)
Pm = (Pm - Pm.mean(0)) / Pm.std(0)
R = (Pm.T @ Pm) / len(Pm)
np.fill_diagonal(R, -9)
fam = K.loc[kin, "manning_family"].values
sets = {}; meta = []
for i, k in enumerate(kin):
    nei = [kin[j] for j in np.argsort(-R[i])[:NNEI]]
    fm = [kin[j] for j in range(len(kin)) if fam[j] == fam[i] and j != i]
    conf = sorted(set(nei) | set(fm))
    ci = [kin.index(c) for c in conf]
    pk = np.nan_to_num(Pv[:, i], nan=0)
    pc = np.nan_to_num(Pv[:, ci], nan=0).max(1) if ci else np.zeros(len(pk))
    full = np.where(pk >= TOP)[0]
    spec = np.where((pk >= TOP) & (pc < MARG))[0]
    cset = np.where((pc >= TOP) & (pk < MARG))[0]
    sets[k] = dict(full=full, spec=spec, conf=cset)
    meta.append(dict(kinase=k, family=fam[i], n_family=len(fm), neighbours=";".join(nei),
                     neighbour_r=";".join(f"{R[i, kin.index(n)]:.2f}" for n in nei), n_conf=len(conf),
                     n_full=len(full), n_spec=len(spec), n_confsites=len(cset), frac_spec=len(spec) / max(len(full), 1)))
meta = pd.DataFrame(meta).set_index("kinase")
meta["dark"] = K.loc[meta.index, "idg_dark_kinase"].values
meta.to_csv(OUT / f"tables/substrate_sets{TAG}.tsv", sep="\t")
np.save(CA / f"sets{TAG}.npy", sets, allow_pickle=True)


def score(idx_lists):
    n_s = Zf.shape[0]
    rows, cols = [], []
    for j, ix in enumerate(idx_lists):
        rows += list(ix); cols += [j] * len(ix)
    I = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(n_s, len(idx_lists)))
    s = (I.T @ Zf).T if False else np.asarray((I.T @ Zf))  # sets x conds
    n = np.asarray(I.T @ obs.astype(np.float32))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = s / np.sqrt(n)
    out[n < NMIN] = np.nan
    return out, n


for typ in ["full", "spec", "conf"]:
    sc, n = score([sets[k][typ] for k in kin])
    pd.DataFrame(sc, index=kin, columns=X.columns).to_parquet(CA / f"score_{typ}{TAG}.parquet")
    pd.DataFrame(n, index=kin, columns=X.columns).to_parquet(CA / f"nsub_{typ}{TAG}.parquet")
print(meta.groupby("dark")[["n_full", "n_spec", "n_confsites", "frac_spec"]].median())
print(meta[meta.dark].sort_values("n_spec").head(20))
