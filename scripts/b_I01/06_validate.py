"""Orthogonal validation of perturbation-derived kinase-kinase associations.
Pairs (q, j): q = query kinase with >=20 specific sites; j any other kinase NOT in q's Manning family and NOT one of
q's 5 motif neighbours. Labels: DepMap co-dependency (top 2% r of all kinase pairs), GTEx co-expression (top 5%),
PPI experimental screens, PPI literature, HPA co-localisation (Jaccard>=0.5), COMPARTMENTS cosine (top 10%).
AUROC per method; null = query-row permutation (scores of query q assigned to another query q'), 300 perms.
Also: 'family-label' baseline = mean label of q's family members with j (guilt-by-family in the orthogonal data)."""
import numpy as np, pandas as pd, sys, os
from pathlib import Path
from scipy.stats import rankdata
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
TAG = os.environ.get("TAG", "")
SET_TAG = os.environ.get("SET_TAG", TAG)
rng = np.random.default_rng(1)
NPERM = 300
M = pd.read_csv(OUT / f"tables/substrate_sets{SET_TAG}.tsv", sep="\t", index_col=0)
T = pd.read_csv(OUT / f"tables/specificity_per_kinase{TAG}.tsv", sep="\t", index_col=0)
kin = list(M.index); nk = len(kin); ki = {k: i for i, k in enumerate(kin)}
A = {nm: pd.read_parquet(CA / f"assoc_{nm}{TAG}.parquet").loc[kin, kin].values for nm in ["resid", "full", "conf", "corev3"]}
labs_raw = {nm: pd.read_parquet(CA / f"lab_{nm}.parquet").reindex(index=kin, columns=kin).values.astype(float)
            for nm in os.environ.get("LABELS", "depmap,depmap_noprox,gtex,ppi_exp,ppi_lit,hpa,compartments").split(",")}
iu = np.triu_indices(nk, 1)


def binarize(nm, X):
    if nm in ("ppi_exp", "ppi_lit"): return X
    v = X[iu]; v = v[~np.isnan(v)]
    thr = {"depmap": np.quantile(v, 0.98), "depmap_noprox": np.quantile(v, 0.98), "gtex": np.quantile(v, 0.95), "hpa": 0.5, "compartments": np.quantile(v, 0.90)}[nm]
    return np.where(np.isnan(X), np.nan, (X >= thr).astype(float))


labs = {nm: binarize(nm, X) for nm, X in labs_raw.items()}
fam = M.family.values
excl = np.zeros((nk, nk), bool)
for i, k in enumerate(kin):
    excl[i, i] = True
    excl[i, fam == fam[i]] = True
    for n in M.loc[k, "neighbours"].split(";"): excl[i, ki[n]] = True
# family-label baseline
famlab = {}
for nm, L in labs.items():
    F = np.full((nk, nk), np.nan)
    for i in range(nk):
        mem = np.where((fam == fam[i]) & (np.arange(nk) != i))[0]
        if len(mem): F[i] = np.nanmean(L[mem], 0)
    famlab[nm] = F
# partial (resid beyond full & conf) association
Ar, Af, Ac = A["resid"], A["full"], A["conf"]
P = np.full_like(Ar, np.nan)
for i in range(nk):
    y = Ar[i]; X = np.c_[np.ones(nk), Af[i], Ac[i]]; m = ~np.isnan(y) & ~np.isnan(X).any(1) & ~excl[i]
    if m.sum() > 20:
        b = np.linalg.lstsq(X[m], y[m], rcond=None)[0]; P[i] = y - X @ b
A["resid_beyond_full_conf"] = P


def auc(s, y):
    m = ~np.isnan(s) & ~np.isnan(y)
    s, y = s[m], y[m]; npos = y.sum(); nneg = len(y) - npos
    if npos < 5: return np.nan, int(npos)
    r = rankdata(s)
    return (r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg), int(npos)


from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
relfdr = np.ones(nk)
_m = T.dark.reindex(kin).fillna(False).values.astype(bool) & T.rel_z.reindex(kin).notna().values
relfdr[_m] = multipletests(norm.sf(T.rel_z.reindex(kin).values[_m]), method="fdr_bh")[1]
hasq = T.n_spec.reindex(kin).fillna(0).values >= 20
dark = T.dark.reindex(kin).fillna(False).values.astype(bool)
qsets = {"dark": np.where(hasq & dark)[0], "nondark": np.where(hasq & ~dark)[0],
         "dark_reliable": np.where(hasq & dark & (relfdr <= 0.05))[0]}
rows = []
for qs, Q in qsets.items():
    for ln, L in labs.items():
        Y = L[Q]; EX = excl[Q]
        y = np.where(EX, np.nan, Y)
        meths = dict(A); meths["family_label"] = famlab[ln]
        for mn, S in meths.items():
            s = np.where(EX, np.nan, S[Q])
            a, npos = auc(s.ravel(), y.ravel())
            if np.isnan(a): continue
            null = []
            if mn != "family_label":
                for _ in range(NPERM):
                    sp_ = np.where(EX, np.nan, S[rng.permutation(Q)])
                    null.append(auc(sp_.ravel(), y.ravel())[0])
                null = np.array(null)
            d = dict(queries=qs, n_queries=len(Q), label=ln, method=mn, auroc=a, n_pos=npos,
                     n_pairs=int((~np.isnan(s) & ~np.isnan(y)).sum()))
            if len(null):
                d.update(null_mean=null.mean(), null_sd=null.std(), p_perm=(1 + (null >= a).sum()) / (1 + NPERM))
            rows.append(d)
        print(qs, ln, flush=True)
R = pd.DataFrame(rows)
R.to_csv(OUT / f"tables/validation_auroc{TAG}{os.environ.get('LTAG', '')}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(R.pivot_table(index=["queries", "label"], columns="method", values="auroc").round(3))
print(R.pivot_table(index=["queries", "label"], columns="method", values="p_perm").round(3))
