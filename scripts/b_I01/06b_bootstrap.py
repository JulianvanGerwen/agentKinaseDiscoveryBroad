"""Bootstrap (over query kinases) CI for AUROC differences resid - full / conf / family_label, reliable dark kinases."""
import numpy as np, pandas as pd, os
from pathlib import Path
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"; TAG = os.environ.get("TAG", "")
M = pd.read_csv(OUT / "tables/substrate_sets.tsv", sep="\t", index_col=0); kin = list(M.index); nk = len(kin); ki = {k: i for i, k in enumerate(kin)}
T = pd.read_csv(OUT / f"tables/specificity_per_kinase{TAG}.tsv", sep="\t", index_col=0).reindex(kin)
A = {n: pd.read_parquet(CA / f"assoc_{n}{TAG}.parquet").loc[kin, kin].values for n in ["resid", "full", "conf"]}
fam = M.family.values
excl = np.zeros((nk, nk), bool)
for i, k in enumerate(kin):
    excl[i, fam == fam[i]] = True; excl[i, i] = True
    for n in M.loc[k, "neighbours"].split(";"): excl[i, ki[n]] = True
m = T.dark.fillna(False).values.astype(bool) & T.rel_z.notna().values & (T.n_spec.fillna(0).values >= 20)
fdr = np.ones(nk); fdr[m] = multipletests(norm.sf(T.rel_z.values[m]), method="fdr_bh")[1]
Q = np.where(m & (fdr <= 0.05))[0]
rng = np.random.default_rng(3); rows = []
for ln in ["depmap_noprox", "compartments"]:
    X = pd.read_parquet(CA / f"lab_{ln}.parquet").reindex(index=kin, columns=kin).values
    v = X[np.triu_indices(nk, 1)]; v = v[~np.isnan(v)]
    L = np.where(np.isnan(X), np.nan, (X >= np.quantile(v, 0.98 if "depmap" in ln else 0.90)).astype(float))
    F = np.full((nk, nk), np.nan)
    for i in range(nk):
        mem = np.where((fam == fam[i]) & (np.arange(nk) != i))[0]
        if len(mem): F[i] = np.nanmean(L[mem], 0)
    S = dict(A); S["family_label"] = F
    def auc(q, s):
        y = np.where(excl[q], np.nan, L[q]).ravel(); x = np.where(excl[q], np.nan, s[q]).ravel()
        k = ~np.isnan(x) & ~np.isnan(y); x, y = x[k], y[k]; p = y.sum()
        r = rankdata(x); return (r[y == 1].sum() - p * (p + 1) / 2) / (p * (len(y) - p))
    obs = {n: auc(Q, s) for n, s in S.items()}
    boots = {n: [] for n in S}
    for _ in range(500):
        qb = rng.choice(Q, len(Q))
        for n, s in S.items(): boots[n].append(auc(qb, s))
    for n in S:
        d = np.array(boots["resid"]) - np.array(boots[n])
        rows.append(dict(label=ln, method=n, n_queries=len(Q), auroc=obs[n], ci_lo=np.percentile(boots[n], 2.5), ci_hi=np.percentile(boots[n], 97.5),
                         d_resid_minus=obs["resid"] - obs[n], d_ci_lo=np.percentile(d, 2.5), d_ci_hi=np.percentile(d, 97.5)))
R = pd.DataFrame(rows); R.to_csv(OUT / f"tables/validation_bootstrap{TAG}.tsv", sep="\t", index=False); print(R.round(3).to_string())
print("queries:", [kin[i] for i in Q])
