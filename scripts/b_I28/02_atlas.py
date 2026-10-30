"""I28 step 2: turnover atlas, matched-random null, split-half reliability (unit-level and study-level), specificity.
Out: tables/turnover_atlas.tsv, tables/group_summary.tsv, tables/reliability.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
sys.path.insert(0, os.path.dirname(__file__)); from lib28 import *
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
rng = np.random.default_rng(1)
L = pd.read_parquet(f"{C28}/L_units.parquet")
S = pd.read_csv(f"{C28}/substrates.tsv", sep="\t")
A = classify_units(L[L.set == "substrate"]).merge(S, on=["group", "site"], how="left")
AR = classify_units(L[L.set == "random"])
SA = load.sites().set_index("site")
A["gene"] = A.site.map(SA.gene); A["window"] = A.site.map(SA.sequence_window); A["psp_regulatory"] = A.site.map(SA.psp_regulatory)
A["n_conditions_compendium"] = A.site.map(SA.n_conditions)

# ---- generic lability: mean z of the site across verified conditions of OTHER groups (excluding studies used for this group)
Z = pd.read_parquet(f"{C26}/Z_lof.parquet").set_index("site")
CU = pd.read_csv(f"{C28}/conds_used.tsv", sep="\t", dtype={"study": str})
lab = {}
for g in A.group.unique():
    st_g = set(L[L.group == g].study)
    oc = CU[CU.worked & (CU.group != g) & ~CU.study.isin(st_g)].conditionID.unique()
    if g == "ERK":
        oc = CU[CU.worked & ~CU.group.isin(["ERK", "MEK", "RAF"]) & ~CU.study.isin(st_g)].conditionID.unique()
    ss = A[A.group == g].site
    lab[g] = Z.reindex(ss)[list(oc)].mean(axis=1)
A["lability_other"] = [lab[g].get(s, np.nan) for g, s in zip(A.group, A.site)]
A.to_csv(f"{R28}/tables/turnover_atlas.tsv", sep="\t", index=False)

# ---- group summary vs matched random
rows = []
for g in sorted(A.group.unique()):
    a = A[(A.group == g) & (A.n_units >= 2)]; r = AR[(AR.group == g) & (AR.n_units >= 2)]
    if len(a) < 5:
        continue
    nst = L[L.group == g].study.nunique()
    def bs(x, f, B=500):
        x = np.asarray(x); v = [f(rng.choice(x, len(x))) for _ in range(B)]; return np.percentile(v, [2.5, 97.5])
    ff = (a.cls == "fast").mean(); fr = (r.cls == "fast").mean()
    lo, hi = bs((a.cls == "fast").astype(float), np.mean)
    rows.append(dict(group=g, n_studies=nst, n_units=L[L.group == g].unit.nunique(), n_sites=len(a),
                     n_psp=(a.source != "pred").sum(), n_pred=(a.source == "pred").sum(),
                     frac_fast=ff, frac_fast_lo=lo, frac_fast_hi=hi, frac_stable=(a.cls == "stable").mean(),
                     frac_fast_random=fr, frac_stable_random=(r.cls == "stable").mean(),
                     excess_fast=ff - fr, median_R=a.R.median(), median_R_random=r.R.median(),
                     frac_fast_psp=(a[a.source != "pred"].cls == "fast").mean(),
                     frac_fast_pred=(a[a.source == "pred"].cls == "fast").mean() if (a.source == "pred").any() else np.nan,
                     frac_fast_gold=(a[a.gold == True].cls == "fast").mean() if (a.gold == True).any() else np.nan,
                     MWU_p=mannwhitneyu(a.R, r.R).pvalue))
GS = pd.DataFrame(rows); GS.to_csv(f"{R28}/tables/group_summary.tsv", sep="\t", index=False)
print(GS.round(3).to_string())


# ---- split-half reliability
def split_half(Lx, level="unit", reps=200, min_n=2):
    """Per rep: for each (group,site) with >=min_n levels, random half split; Spearman/Pearson of half means, after
    removing group means (within-family reliability). Spearman-Brown corrected."""
    key = "unit" if level == "unit" else "study"
    if level == "study":
        Lx = Lx.groupby(["group", "site", "study"], observed=True).z.mean().reset_index()
    cnt = Lx.groupby(["group", "site"]).size()
    keep = cnt[cnt >= min_n].index
    Lx = Lx.set_index(["group", "site"]).loc[keep].reset_index()
    Lx["k"] = Lx.groupby(["group", "site"]).cumcount(); Lx["n"] = Lx.groupby(["group", "site"]).z.transform("size")
    res = []
    for r in range(reps):
        u = rng.random(len(Lx)); Lx["rk"] = Lx.assign(u=u).groupby(["group", "site"]).u.rank()
        h = np.where(Lx.rk <= Lx.n / 2, "a", "b")
        P = Lx.assign(h=h).groupby(["group", "site", "h"]).z.mean().unstack()
        P = P.dropna()
        P = P - P.groupby(level=0).transform("mean")
        rho = spearmanr(P.a, P.b)[0]; pr = pearsonr(P.a, P.b)[0]
        res.append((rho, pr))
    res = np.array(res)
    rho = res[:, 0].mean(); pr = res[:, 1].mean()
    return dict(n_sites=len(keep), rho=rho, rho_lo=np.percentile(res[:, 0], 2.5), rho_hi=np.percentile(res[:, 0], 97.5),
                pearson=pr, spearman_brown=2 * pr / (1 + pr))


rel = []
for setn in ["substrate", "random"]:
    for level in ["unit", "study"]:
        for g in [None] + sorted(A.group.unique()):
            Lx = L[(L.set == setn)] if g is None else L[(L.set == setn) & (L.group == g)]
            try:
                d = split_half(Lx[["group", "site", "unit", "study", "z"]], level, reps=100 if g else 200)
            except Exception:
                continue
            if d["n_sites"] < 15:
                continue
            rel.append(dict(set=setn, level=level, group=g or "ALL", **d))
# cross-study reliability with the generic-lability covariate partialled out (substrates, study level)
Ls = L[L.set == "substrate"].merge(A[["group", "site", "lability_other"]], on=["group", "site"])
Ls = Ls[Ls.lability_other.notna()]
Ls["z"] = Ls.z - np.polyval(np.polyfit(Ls.lability_other, Ls.z, 1), Ls.lability_other)
d = split_half(Ls[["group", "site", "unit", "study", "z"]], "study", reps=200)
rel.append(dict(set="substrate_lability_adjusted", level="study", group="ALL", **d))
RL = pd.DataFrame(rel); RL.to_csv(f"{R28}/tables/reliability.tsv", sep="\t", index=False)
print(RL.round(3).to_string())

# specificity: is R correlated with lability in other-group inhibitions?
x = A.dropna(subset=["lability_other"]); x = x[x.n_units >= 2]
print("R vs lability_other spearman", spearmanr(x.R, x.lability_other))
print(A[A.n_units >= 2].groupby("group").cls.value_counts().unstack(fill_value=0).to_string())
