"""Sensitivity: DepMap validation restricted to kinases with non-trivial Chronos variation (SD above the kinase
33rd percentile), since co-dependency among non-expressed (tissue-restricted) genes may be noise-driven."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests
E = Path("<HACKATHON_ROOT>/kinaseSubMain/external"); OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
M = pd.read_csv(OUT / "tables/substrate_sets.tsv", sep="\t", index_col=0); kin = list(M.index); nk = len(kin); ki = {k: i for i, k in enumerate(kin)}
D = pd.read_csv(E / "depmap/DepMap_26Q1_chronos_gene_effect.csv", index_col=0)
D.columns = [c.split(" (")[0] for c in D.columns]
sd = D.reindex(columns=kin).std().reindex(kin)
sd.to_csv(CA / "chronos_sd.tsv", sep="\t")
rng = np.random.default_rng(7)
T = pd.read_csv(OUT / "tables/specificity_per_kinase.tsv", sep="\t", index_col=0).reindex(kin)
A = {n: pd.read_parquet(CA / f"assoc_{n}.parquet").loc[kin, kin].values for n in ["resid", "full", "conf"]}
fam = M.family.values; excl = np.zeros((nk, nk), bool)
for i, k in enumerate(kin):
    excl[i, fam == fam[i]] = True; excl[i, i] = True
    for n in M.loc[k, "neighbours"].split(";"): excl[i, ki[n]] = True
m = T.dark.fillna(False).values.astype(bool) & T.rel_z.notna().values & (T.n_spec.fillna(0).values >= 20)
fdr = np.ones(nk); fdr[m] = multipletests(norm.sf(T.rel_z.values[m]), method="fdr_bh")[1]
X = pd.read_parquet(CA / "lab_depmap_noprox.parquet").reindex(index=kin, columns=kin).values
v = X[np.triu_indices(nk, 1)]; v = v[~np.isnan(v)]; L = np.where(np.isnan(X), np.nan, (X >= np.quantile(v, 0.98)).astype(float))
def auc(Q, s, ok, Qs=None):
    Qs = Q if Qs is None else Qs
    y = np.where(excl[Q] | ~ok[Q], np.nan, L[Q]).ravel(); x = np.where(excl[Q], np.nan, s[Qs]).ravel()
    k = ~np.isnan(x) & ~np.isnan(y); x, y = x[k], y[k]; p = y.sum()
    if p < 5: return np.nan, 0
    r = rankdata(x); return (r[y == 1].sum() - p * (p + 1) / 2) / (p * (len(y) - p)), int(p)
rows = []
for q in [0.0, 0.33, 0.5]:
    thr = sd.quantile(q); e = (sd >= thr).values
    ok = e[:, None] & e[None, :]
    Q = np.where(m & (fdr <= 0.05) & e)[0]
    for n, s in A.items():
        a, p = auc(Q, s, ok)
        null = [auc(Q, s, ok, rng.permutation(Q))[0] for _ in range(200)] if n == "resid" else []
        rows.append(dict(sd_quantile=q, method=n, n_queries=len(Q), n_pos=p, auroc=a, p_perm=(1 + np.sum(np.array(null) >= a)) / 201 if null else np.nan))
        print(f"sd>=q{q}", n, "nq", len(Q), "npos", p, round(a, 3), ("p_perm=%.3f" % ((1 + np.sum(np.array(null) >= a)) / 201)) if null else "")
pd.DataFrame(rows).to_csv(OUT / "tables/validation_depmap_expressed.tsv", sep="\t", index=False)
