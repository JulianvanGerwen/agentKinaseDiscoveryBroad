"""Audit H2b: pThr effect without the McCloy 2015 mitotic-exit study (26055452); per-study T vs S effects."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import mannwhitneyu
sys.path.insert(0, "/workspace/scripts/b_I28"); from lib28 import classify_units, C28
OUT = "/workspace/results/audit_H"
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]
L = pd.read_parquet(f"{C28}/L_units.parquet"); L = L[L.group.isin(MAIN)]
F = pd.read_parquet(f"{C28}/features.parquet")
L["res"] = L.site.str.split("_").str[1].str[0]; L = L[L.res.isin(["S", "T"])]
rows = []
for (setn, g, st), d in L.groupby(["set", "group", "study"]):
    x = d.groupby("site").agg(z=("z", "mean"), res=("res", "first"))
    t, s = x.z[x.res == "T"], x.z[x.res == "S"]
    if len(t) >= 5 and len(s) >= 5:
        rows.append(dict(set=setn, group=g, study=st, n_T=len(t), n_S=len(s), med_T=t.median(), med_S=s.median(),
                         frac_dec_T=(t <= -1).mean(), frac_dec_S=(s <= -1).mean(), p=mannwhitneyu(t, s).pvalue))
P = pd.DataFrame(rows); P.to_csv(f"{OUT}/h2_per_study_T_vs_S.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(P.round(3).to_string())
out = []
for drop in [None, "26055452"]:
    for setn in ["substrate", "random"]:
        A = classify_units(L[(L.set == setn) & (L.study != drop)]); A = A[A.n_units >= 2].join(F, on="site")
        d = A[A.cls.isin(["fast", "stable"])].dropna(subset=["P_plus1", "basic_minus3", "n_acidic_pm5"])
        d = d.assign(isT=d.is_T, prot=d.site.str.split("_").str[0])
        X = sm.add_constant(pd.concat([d[["isT", "P_plus1", "basic_minus3", "n_acidic_pm5"]], np.log2(d.n_units).rename("lnu"),
                                       pd.get_dummies(d.group, drop_first=True, dtype=float)], axis=1))
        m = sm.Logit((d.cls == "fast").astype(float), X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": d.prot.astype("category").cat.codes})
        out.append(dict(dropped=drop or "none", set=setn, n=len(d), n_T=int(d.isT.sum()), OR_T=np.exp(m.params.isT),
                        lo=np.exp(m.params.isT - 1.96 * m.bse.isT), hi=np.exp(m.params.isT + 1.96 * m.bse.isT), p=m.pvalues.isT,
                        OR_acidic=np.exp(m.params.n_acidic_pm5), p_acidic=m.pvalues.n_acidic_pm5,
                        OR_basic3=np.exp(m.params.basic_minus3), p_basic3=m.pvalues.basic_minus3))
O = pd.DataFrame(out); O.to_csv(f"{OUT}/h2_drop_mccloy.tsv", sep="\t", index=False); print(O.round(3).to_string())
