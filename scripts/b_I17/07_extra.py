"""Follow-ups: (1) most-potent-site (top decile pEC50, down) vs pKd_target with drug-label shuffle;
(2) target-proximal (self + curated substrate) shift per experiment; (3) shift model adding residue type and
per-drug pred_direct effects (leave-one-drug-out); (4) binding mode (Klaeger S1) vs shift; (5) reproducibly steep sites.
Outputs: results/b_I17/extra_*.tsv
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.formula.api as smf
OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"
A = pd.read_csv(f"{RES}/shift_sites.tsv", sep="\t")
F = pd.read_parquet(f"{RES}/curves_all.parquet")
X = pd.read_csv(f"{RES}/shift_by_experiment.tsv", sep="\t")
s1 = pd.read_excel("/workspace/data/Klaeger/aan4368_Table_S1.xlsx", sheet_name="Inhibitor annotation").set_index("Drug")
rng = np.random.default_rng(7)
A["residue"] = A.site.str.extract(r"_([STY])\d+$")[0]
A["direct"] = A.ks_direct | A.self
dn = A[A.direction == "down"].copy()
X = X[X.n_down >= 10].copy()

# (1) top-decile potency
top = dn.groupby("exp_id").pEC50.quantile(0.9).rename("p90_pEC50_down")
X = X.join(top, on="exp_id")
prox = dn[dn.direct].groupby("exp_id").pEC50.median().rename("median_pEC50_direct")
X = X.join(prox, on="exp_id")
rows = []
for col in ["median_pEC50_down", "p90_pEC50_down", "median_pEC50_direct"]:
    d = X.dropna(subset=[col])
    dd = d.groupby("agent").agg(pKd=("pKd_target", "first"), y=(col, "median"))
    r = spearmanr(dd.pKd, dd.y)[0]
    null = np.array([spearmanr(rng.permutation(dd.pKd.values), dd.y)[0] for _ in range(5000)])
    rows.append(dict(readout=col, level="drug", n=len(dd), rho=r, p_shuffle=(np.sum(null >= r) + 1) / 5001,
                     median_shift=(dd.y - dd.pKd).median()))
R1 = pd.DataFrame(rows); R1.to_csv(f"{RES}/extra_potency_vs_pKd.tsv", sep="\t", index=False); print(R1.round(3).to_string())

# (4) binding mode / type
X["binding_mode"] = X.agent.map(lambda a: None)
km = {"afatinib": "Afatinib", "dasatinib": "Dasatinib", "gefitinib": "Gefitinib", "imatinib": "Imatinib", "MK2206": "MK-2206",
      "PD0325901": "PD-325901", "nintedanib": "Nintedanib", "refametinib": "Refametinib", "selumetinib": "Selumetinib",
      "lapatinib": "Lapatinib", "CUDC101": "CUDC-101", "AZD4547": "AZD-4547"}
X["binding_mode"] = [s1.loc[km[a], "Binding mode"] if a in km else None for a in X.agent]
X["binding_type"] = [s1.loc[km[a], "Binding type"] if a in km else None for a in X.agent]
X["allosteric"] = X.agent.isin(["MK2206", "PD0325901", "refametinib", "selumetinib"])
X.to_csv(f"{RES}/extra_shift_by_experiment.tsv", sep="\t", index=False)
dd = X.groupby("agent").agg(shift=("median_shift_down", "median"), allo=("allosteric", "first"), mode=("binding_mode", "first"))
u = mannwhitneyu(dd["shift"][dd.allo.astype(bool)], dd["shift"][~dd.allo.astype(bool)], alternative="greater")
print(dd.sort_values("shift").round(2).to_string()); print("allosteric vs ATP-site drug-level MWU p", u.pvalue)

# (3) model with residue + leave-one-drug-out for pred_direct
dn["dist_c"] = dn.string_dist.clip(upper=4)
form = "shift ~ direct + pred_direct + C(dist_c) + C(residue) + pmembrane + nucleus + C(exp_id)"
m = smf.ols(form, data=dn.dropna(subset=["dist_c", "residue"])).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dn.dropna(subset=["dist_c", "residue"]).agent)[0]})
print(m.summary().tables[1])
lodo = []
for a in dn.agent.unique():
    d = dn[(dn.agent != a)].dropna(subset=["dist_c", "residue"])
    mm = smf.ols(form, data=d).fit()
    lodo.append(dict(left_out=a, coef_pred_direct=mm.params.get("pred_direct[T.True]"), coef_direct=mm.params.get("direct[T.True]")))
per = []
for ex, d in dn.groupby("exp_id"):
    if d.pred_direct.sum() >= 3 and (~d.pred_direct).sum() >= 3:
        per.append(dict(exp_id=ex, agent=d.agent.iloc[0], n_pred=int(d.pred_direct.sum()),
                        diff=d.loc[d.pred_direct, "shift"].median() - d.loc[~d.pred_direct, "shift"].median(),
                        frac_Y_pred=(d.loc[d.pred_direct, "residue"] == "Y").mean()))
P = pd.DataFrame(per); L = pd.DataFrame(lodo)
P.to_csv(f"{RES}/extra_pred_direct_per_experiment.tsv", sep="\t", index=False)
L.to_csv(f"{RES}/extra_pred_direct_lodo.tsv", sep="\t", index=False)
pd.DataFrame(dict(coef=m.params, se=m.bse, p=m.pvalues)).filter(regex="^(?!C\\(exp)", axis=0).to_csv(f"{RES}/extra_shift_model_residue.tsv", sep="\t")
print(P.round(2).to_string()); print(L.round(3).to_string())
print("experiments with pred_direct diff>0:", (P["diff"] > 0).sum(), "/", len(P))

# (5) reproducibly steep / shallow sites (dasatinib A549 replicate pair) and up vs down Hill
a = F[(F.exp_id == "DS0015_Dasatinib0nM_0") & F.sig].set_index("site"); b = F[(F.exp_id == "DS0015_Dasatinib0nM_1") & F.sig].set_index("site")
c = a.index.intersection(b.index)
st = pd.DataFrame(dict(hill_a=a.loc[c, "hill"], hill_b=b.loc[c, "hill"], pEC50_a=a.loc[c, "pEC50"], pEC50_b=b.loc[c, "pEC50"], delta_a=a.loc[c, "delta"]))
st["steep_both"] = (st.hill_a >= 2) & (st.hill_b >= 2)
S = pd.read_parquet(f"{OUT}/sites_min.parquet").set_index("site")
st["gene"] = S.gene.reindex(st.index).values
st.to_csv(f"{RES}/extra_dasatinib_A549_replicate_hill.tsv", sep="\t")
print("both-sig sites", len(st), "steep in both", st.steep_both.sum(), "expected if indep", ((st.hill_a >= 2).mean() * (st.hill_b >= 2).mean() * len(st)).round(1))
print(st[st.steep_both].sort_values("pEC50_a").round(2).to_string())
hu = A.groupby("direction").hill.describe(); print(hu)
