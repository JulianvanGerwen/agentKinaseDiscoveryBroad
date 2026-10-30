"""Audit H2: is the pThr->fast association a composition artefact?
Analyses on S/T substrates of the 5 main families (n_units>=2):
 (a) logistic fast vs stable: group FE (I28) -> + motif covariates (P+1, basic-3, acidic, motif pct, source) -> + study-coverage pattern
 (b) continuous -R, within (group x protein) fixed effects (demeaned OLS), cluster bootstrap by protein
 (c) paired within-protein: proteins with >=1 pS and >=1 pT substrate in same family: mean R(T) - mean R(S), sign test
 (d) same for random matched sites (generic residue effect under inhibition?) and per family / per lab
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import binomtest, wilcoxon, mannwhitneyu
sys.path.insert(0, "/workspace/scripts/b_I28"); from lib28 import classify_units, C28
OUT = "/workspace/results/audit_H"; rng = np.random.default_rng(11)
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]
A = pd.read_csv("/workspace/results/b_I28/tables/turnover_atlas.tsv", sep="\t")
A = A[A.group.isin(MAIN) & (A.n_units >= 2)].copy()
F = pd.read_parquet(f"{C28}/features.parquet")
L = pd.read_parquet(f"{C28}/L_units.parquet"); L = L[L.group.isin(MAIN)]
AR = classify_units(L[L.set == "random"]); AR = AR[AR.n_units >= 2].copy()
LAB = {"31959955": "Cutillas", "26060313": "Cutillas", "36926954": "Kuster", "27794612": "Kuster", "28209619": "Kuster",
       "30901187": "Kuster", "25850435": "Ahn", "29255136": "Ahn", "30459283": "Olsen", "32034161": "Olsen"}
rows = []


def prep(D):
    D = D.join(F.drop(columns=["psp_regulatory"]), on="site"); D["res"] = D.site.str.split("_").str[1].str[0]
    D = D[D.res.isin(["S", "T"])].copy(); D["prot"] = D.site.str.split("_").str[0]
    D["isT"] = (D.res == "T").astype(float)
    # which studies measured the site (coverage pattern) -> set of labs
    labs = L.assign(lab=L.study.map(lambda s: LAB.get(s, "L" + s))).groupby(["group", "site"]).lab.apply(lambda x: "|".join(sorted(set(x))))
    D["labset"] = labs.reindex(pd.MultiIndex.from_frame(D[["group", "site"]])).values
    D["lognu"] = np.log2(D.n_units)
    return D


S = prep(A); R = prep(AR)
print("pT fraction by group (substrates):", S.groupby("group").isT.mean().round(3).to_dict())
print("pT fraction by group (random):", R.groupby("group").isT.mean().round(3).to_dict())
print("pT fraction among P+1 vs not (substrates):", S.groupby("P_plus1").isT.mean().round(3).to_dict())


def logit(D, covs, fe, tag):
    d = D[D.cls.isin(["fast", "stable"])].dropna(subset=covs)
    X = pd.concat([d[["isT"] + covs]] + [pd.get_dummies(d[f], drop_first=True, dtype=float, prefix=f) for f in fe], axis=1)
    X = X.loc[:, X.std() > 0]; X = sm.add_constant(X)
    y = (d.cls == "fast").astype(float)
    m = sm.Logit(y, X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": d.prot.astype("category").cat.codes}, maxiter=200)
    b, se = m.params.isT, m.bse.isT
    rows.append(dict(analysis=tag, n=len(d), n_T=int(d.isT.sum()), OR_T=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues.isT))


motif = ["P_plus1", "basic_minus3", "n_acidic_pm5", "n_basic_pm5"]
for D, nm in [(S, "substrate"), (R, "random")]:
    logit(D, [], ["group"], f"{nm}: group FE (I28-like, no coverage)")
    logit(D, ["lognu"], ["group"], f"{nm}: group FE + n_units")
    logit(D, ["lognu"] + motif, ["group"], f"{nm}: + motif covariates")
    logit(D, ["lognu"] + motif + ["idr", "plddt"], ["group", "labset"], f"{nm}: + motif + IDR/pLDDT + lab-coverage-pattern FE")
    for g in MAIN:
        try:
            logit(D[D.group == g], ["lognu"] + motif, [], f"{nm}: {g} only, + motif")
        except Exception as e:
            print(g, e)
    # P+1 strata
    for p1 in [0, 1]:
        try: logit(D[D.P_plus1 == p1], ["lognu", "basic_minus3", "n_acidic_pm5"], ["group"], f"{nm}: P+1={p1} stratum")
        except Exception as e: print(e)
    # family-specific: CDK, ERK (proline-directed) the T effect
# continuous: within group x protein FE
def within_prot(D, tag, B=500):
    d = D.copy(); d["gp"] = d.group + "|" + d.prot
    k = d.groupby("gp").isT.agg(["min", "max", "size"]); inf = k[(k["min"] == 0) & (k["max"] == 1)].index
    d = d[d.gp.isin(inf)].copy()
    y = -d.R.clip(-6, 4); cov = d[["isT", "P_plus1", "basic_minus3", "n_acidic_pm5", "lognu"]].astype(float)
    yd = y - y.groupby(d.gp).transform("mean"); Xd = cov - cov.groupby(d.gp).transform("mean")
    def fit(ix):
        return np.linalg.lstsq(Xd.values[ix], yd.values[ix], rcond=None)[0][0]
    b = fit(np.arange(len(d)))
    gps = d.gp.unique(); pos = {g: np.where(d.gp.values == g)[0] for g in gps}
    bs = [fit(np.concatenate([pos[g] for g in rng.choice(gps, len(gps))])) for _ in range(B)]
    # paired sign test
    pm = d.groupby(["gp", "res"]).R.mean().unstack()
    diff = (pm["S"] - pm["T"])  # >0 : T more negative (faster)
    rows.append(dict(analysis=f"{tag}: within group x protein, -R ~ isT + motif (OLS FE)", n=len(d), n_T=int(d.isT.sum()),
                     beta_negR_T=b, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), n_proteins=len(gps),
                     frac_prot_T_faster=(diff > 0).mean(), sign_p=binomtest(int((diff > 0).sum()), int((diff != 0).sum())).pvalue))


within_prot(S, "substrate"); within_prot(R, "random")
for g in ["CDK1/2/5/9", "ERK", "AKT"]:
    within_prot(S[S.group == g], f"substrate {g}")
# same-protein +-15 neighbours only
O = pd.DataFrame(rows); O.to_csv(f"{OUT}/h2_pthr.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(O.round(3).to_string())
# per-lab: is T effect present in each lab's data separately?
LS = L[L.set == "substrate"].copy(); LS["lab"] = LS.study.map(lambda s: LAB.get(s, "L" + s)); LS["isT"] = LS.site.str.split("_").str[1].str[0].eq("T")
LS = LS[LS.site.str.split("_").str[1].str[0].isin(["S", "T"])]
x = LS.groupby(["lab", "group", "site"]).agg(z=("z", "mean"), isT=("isT", "first")).reset_index()
x["zc"] = x.z - x.groupby(["lab", "group"]).z.transform("median")
print("per-lab median within-family-centred z, T vs S:")
print(x.groupby(["lab", "isT"]).zc.median().unstack().assign(n=x.groupby("lab").size()).round(3).to_string())
