"""Audit H1: lab structure, leave-one-lab-out, site-noise confounds, unrelated-perturbation reproducibility.
Out: results/audit_H/h1_*.tsv, cache in results/audit_H/cache/
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
sys.path.insert(0, "/workspace/scripts/b_I28"); from lib28 import classify_units, C28, C26
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/audit_H"; os.makedirs(f"{OUT}/cache", exist_ok=True)
rng = np.random.default_rng(7)
LAB = {"31959955": "Cutillas", "26060313": "Cutillas", "36926954": "Kuster", "27794612": "Kuster", "28209619": "Kuster",
       "30901187": "Kuster", "25850435": "Ahn", "29255136": "Ahn", "30459283": "Olsen", "32034161": "Olsen"}
lab_of = lambda s: LAB.get(s, "L" + s)
L = pd.read_parquet(f"{C28}/L_units.parquet")
L["lab"] = L.study.map(lab_of)
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]
L = L[L.group.isin(MAIN)]


def split_half(Lx, key, reps=200, fixed=None):
    """Within-family-centred split-half Spearman of `key`-level means. fixed: function key->'a'/'b' (deterministic split)."""
    Lx = Lx.groupby(["group", "site", key], observed=True).z.mean().reset_index()
    cnt = Lx.groupby(["group", "site"]).size(); keep = cnt[cnt >= 2].index
    Lx = Lx.set_index(["group", "site"]).loc[keep].reset_index()
    if len(keep) < 15:
        return dict(n_sites=len(keep), rho=np.nan, lo=np.nan, hi=np.nan)
    res = []
    for r in range(1 if fixed else reps):
        if fixed:
            h = Lx[key].map(fixed).values
        else:
            u = rng.random(len(Lx)); rk = Lx.assign(u=u).groupby(["group", "site"]).u.rank()
            n = Lx.groupby(["group", "site"]).z.transform("size"); h = np.where(rk <= n / 2, "a", "b")
        P = Lx.assign(h=h).groupby(["group", "site", "h"]).z.mean().unstack().dropna()
        if len(P) < 15: return dict(n_sites=len(P), rho=np.nan, lo=np.nan, hi=np.nan)
        P = P - P.groupby(level=0).transform("mean")
        res.append(spearmanr(P.a, P.b)[0])
    return dict(n_sites=len(P), rho=np.mean(res), lo=np.percentile(res, 2.5), hi=np.percentile(res, 97.5))


rows = []
def add(tag, Lx, key, **kw):
    for setn in ["substrate", "random"]:
        for g in [None] + MAIN:
            d = Lx[(Lx.set == setn) & ((Lx.group == g) if g else True)]
            r = split_half(d, key, reps=100, **kw)
            rows.append(dict(analysis=tag, set=setn, group=g or "ALL", **r))

print(L.groupby(["group", "lab"]).study.nunique().unstack(fill_value=0).to_string())
add("study_split (I28 replicate)", L, "study")
add("lab_split", L, "lab")
add("Cutillas_vs_rest (fixed)", L, "lab", fixed=lambda x: "a" if x == "Cutillas" else "b")
for lab in ["Cutillas", "Kuster", "Ahn", "Olsen"]:
    add(f"study_split_drop_{lab}", L[L.lab != lab], "study")
    add(f"lab_split_drop_{lab}", L[L.lab != lab], "lab")
H = pd.DataFrame(rows); H.to_csv(f"{OUT}/h1_lab_splits.tsv", sep="\t", index=False)
print(H.pivot_table(index=["analysis", "group"], columns="set", values="rho").round(3).to_string())
print(H.pivot_table(index=["analysis", "group"], columns="set", values="n_sites").to_string())

# ---------------- site-noise / dynamism from compendium (non-inhibitor conditions)
Cn = load.conditions()
uncl = ["other", "stress", "disease_vs_control", "genetic", "time_course", "cytokine", "growth_factor"]
sites_all = sorted(set(L.site))
cache = f"{OUT}/cache/Zunrel.parquet"
CU = pd.read_csv(f"{C28}/conds_used.tsv", sep="\t", dtype={"study": str})
Zl = pd.read_parquet(f"{C26}/Z_lof.parquet").set_index("site")
# pool for re-matched randoms: all sites measured in verified units of each group
pool_M = {}
for g in MAIN:
    w = CU[(CU.group == g) & CU.worked]
    M = pd.DataFrame({u: Zl[list(d.conditionID)].mean(axis=1) for u, d in w.groupby("unit")})
    pool_M[g] = M
pool_sites = sorted(set().union(*[set(M.index[M.notna().sum(axis=1) > 0]) for M in pool_M.values()]))
if not os.path.exists(cache):
    cc = Cn[Cn.perturbation_class.isin(uncl) & Cn.qc_flag.eq("ok")]
    X = load.sites_matrix(cc.conditionID, sites=pool_sites)
    med = X.median(axis=0); mad = (X - med).abs().median(axis=0) * 1.4826
    Zu = ((X - med) / mad.clip(lower=0.05)).astype(np.float32)
    Zu.reset_index().to_parquet(cache)
Zu = pd.read_parquet(cache).set_index("site")
cls_of = Cn.set_index("conditionID").perturbation_class.reindex(Zu.columns)
st_of = Cn.set_index("conditionID").study_pmid.reindex(Zu.columns)
# noise: robust SD of z over 'quiet' classes (other/genetic/disease/time_course excluded growth factor & cytokine & stress)
quiet = cls_of.index[cls_of.isin(["other", "genetic", "disease_vs_control", "time_course"])]
Zq = Zu[quiet]
nq = Zq.notna().sum(axis=1)
noise = ((Zq.sub(Zq.median(axis=1), axis=0)).abs().median(axis=1) * 1.4826).where(nq >= 10)
dyn_all = Zu.abs().median(axis=1).where(Zu.notna().sum(axis=1) >= 10)       # generic dynamism (all non-inhibitor)
SA = load.sites().set_index("site")
ncond = np.log10(SA.n_conditions.reindex(Zu.index).astype(float))
SF = pd.DataFrame({"noise_quiet": noise, "dyn_unrel": dyn_all, "log_ncond": ncond, "n_quiet": nq})
SF.to_csv(f"{OUT}/cache/site_noise.tsv", sep="\t")

A = pd.read_csv("/workspace/results/b_I28/tables/turnover_atlas.tsv", sep="\t")
A = A[A.group.isin(MAIN) & (A.n_units >= 2)].join(SF, on="site")
AR = classify_units(L[L.set == "random"]); AR = AR[AR.n_units >= 2].join(SF, on="site")
out2 = []
for nm, d in [("substrate", A), ("random", AR)]:
    for f in ["noise_quiet", "dyn_unrel", "log_ncond"]:
        dd = d.dropna(subset=[f])
        rr = [spearmanr(x[f], -x.R)[0] for _, x in dd.groupby("group") if len(x) > 20]
        out2.append(dict(set=nm, feature=f, mean_within_group_rho_with_negR=np.mean(rr), n=len(dd),
                         median_fast=dd[dd.cls == "fast"][f].median(), median_stable=dd[dd.cls == "stable"][f].median()))
N2 = pd.DataFrame(out2); N2.to_csv(f"{OUT}/h1_noise_vs_R.tsv", sep="\t", index=False); print(N2.round(3).to_string())

# ---------------- re-matched random sites (ST/Y, n units, noise decile, dynamism decile, coverage decile)
Lsub = L[L.set == "substrate"]
SFq = SF.dropna(subset=["noise_quiet", "dyn_unrel"]).copy()
for f in ["noise_quiet", "dyn_unrel", "log_ncond"]:
    SFq[f + "_b"] = pd.qcut(SFq[f].rank(method="first"), 10, labels=False)
S28 = pd.read_csv(f"{C28}/substrates.tsv", sep="\t")
newr = []
for g in MAIN:
    M = pool_M[g]; nu = M.notna().sum(axis=1); nu = nu[nu > 0]
    fam = set(S28[S28.group == g].site)
    pool = SFq.loc[SFq.index.intersection(nu.index.difference(list(fam)))].copy()
    pool["k"] = np.minimum(nu.reindex(pool.index).values, 6); pool["st"] = pool.index.str.split("_").str[1].str[0].isin(["S", "T"])
    key = ["k", "st", "noise_quiet_b", "dyn_unrel_b", "log_ncond_b"]
    grp = pool.groupby(key).groups
    subs = [s for s in Lsub[Lsub.group == g].site.unique() if s in SFq.index and s in nu.index]
    picked = []
    for s in subs:
        kk = (min(nu[s], 6), s.split("_")[1][0] in "ST", *SFq.loc[s, ["noise_quiet_b", "dyn_unrel_b", "log_ncond_b"]].astype(int))
        c = grp.get(tuple(kk))
        if c is not None and len(c):
            picked += list(rng.choice(np.array(c), min(3, len(c)), replace=False))
    picked = sorted(set(picked))
    x = M.loc[picked].stack().dropna().rename("z").reset_index(); x.columns = ["site", "unit", "z"]
    x["group"] = g; x["study"] = x.unit.str.split("|").str[0]; x["set"] = "random_noise_matched"
    newr.append(x)
    print(g, "matched substrates", len(subs), "picked", len(picked))
LR = pd.concat(newr); LR["lab"] = LR.study.map(lab_of)
LS_m = Lsub[Lsub.site.isin(SFq.index)]
rows3 = []
for nm, d in [("substrate (noise-annotated)", LS_m), ("random_noise_matched", LR), ("random_I28", L[(L.set == "random") & L.site.isin(SFq.index)])]:
    for key in ["study", "lab"]:
        for g in [None] + MAIN:
            r = split_half(d[(d.group == g) if g else np.ones(len(d), bool)], key, reps=100)
            rows3.append(dict(set=nm, level=key, group=g or "ALL", **r))
    ARm = classify_units(d); ARm = ARm[ARm.n_units >= 2]
    rows3.append(dict(set=nm, level="frac_fast", group="ALL", n_sites=len(ARm), rho=(ARm.cls == "fast").mean()))
R3 = pd.DataFrame(rows3); R3.to_csv(f"{OUT}/h1_noise_matched_random.tsv", sep="\t", index=False)
print(R3.pivot_table(index=["level", "group"], columns="set", values="rho").round(3).to_string())

# ---------------- substrates: reproducibility of response to UNRELATED perturbations (study split), vs to kinase inhibitors
Zu_s = Zu.reindex(sorted(set(Lsub.site)))
lu = []
for st, cols in pd.Series(Zu.columns, index=Zu.columns).groupby(st_of):
    m = Zu_s[list(cols)].mean(axis=1); a = Zu_s[list(cols)].abs().mean(axis=1)
    lu.append(pd.DataFrame({"site": m.index, "study": st, "z": m.values, "absz": a.values}))
LU = pd.concat(lu).dropna()
fam_sites = Lsub[["group", "site"]].drop_duplicates()
LU = fam_sites.merge(LU, on="site"); LU["lab"] = LU.study.map(lab_of)
# matched random (I28) for same
Lr = L[L.set == "random"]
LUr = Lr[["group", "site"]].drop_duplicates().merge(pd.concat([pd.DataFrame({"site": Zu.index, "study": st,
      "z": Zu[list(c)].mean(axis=1).values, "absz": Zu[list(c)].abs().mean(axis=1).values}) for st, c in
      pd.Series(Zu.columns, index=Zu.columns).groupby(st_of)]).dropna(), on="site")
rows4 = []
for nm, d in [("substrate", LU), ("random_I28", LUr)]:
    for val in ["z", "absz"]:
        dd = d.assign(z=d[val])
        for g in [None] + MAIN:
            r = split_half(dd[(dd.group == g) if g else np.ones(len(dd), bool)], "study", reps=50)
            rows4.append(dict(set=nm, quantity=f"unrelated_perturbations_{val}", group=g or "ALL", **r))
R4 = pd.DataFrame(rows4); R4.to_csv(f"{OUT}/h1_unrelated_reproducibility.tsv", sep="\t", index=False)
print(R4.pivot_table(index=["quantity", "group"], columns="set", values="rho").round(3).to_string())
# does generic unrelated-perturbation signed response predict R?
uz = LU.groupby(["group", "site"]).agg(z_unrel=("z", "mean"), absz_unrel=("absz", "mean")).reset_index()
Ax = A.merge(uz, on=["group", "site"])
print("within-group spearman(-R, z_unrel):", {g: round(spearmanr(-d.R, d.z_unrel)[0], 3) for g, d in Ax.groupby("group")})
print("within-group spearman(-R, absz_unrel):", {g: round(spearmanr(-d.R, d.absz_unrel)[0], 3) for g, d in Ax.groupby("group")})

# ---------------- reliability after residualising z on noise/dynamism/coverage (substrates)
Lz = Lsub.merge(SF, left_on="site", right_index=True).dropna(subset=["noise_quiet", "dyn_unrel", "log_ncond"])
Xd = np.column_stack([np.ones(len(Lz)), Lz.noise_quiet, Lz.dyn_unrel, Lz.log_ncond])
beta = np.linalg.lstsq(Xd, Lz.z.values, rcond=None)[0]; Lz["z"] = Lz.z - Xd @ beta + beta[0]
rr = [dict(analysis="substrate z residualised on noise+dynamism+coverage", level=k, **split_half(Lz, k, reps=100)) for k in ["study", "lab"]]
pd.DataFrame(rr).to_csv(f"{OUT}/h1_resid_noise.tsv", sep="\t", index=False); print(pd.DataFrame(rr).round(3))
