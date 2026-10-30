"""Perturbation contexts for each kinase's specific residual (spec - beta*conf; raw, NOT study-centred, since every
condition is already a contrast vs its own control). Sign flipped for reversed contrasts.
Score standardised per kinase by robust SD over all conditions. Context effect = mean over studies of study-mean;
inference: t-test over study means (>=3 studies), sign consistency; BH over all tests with >=3 studies.
Also computes the same for conf (family/neighbour) and full scores, to show the residual effect is not the family's."""
import numpy as np, pandas as pd, sys, os
from pathlib import Path
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
TAG = os.environ.get("TAG", "")
SET_TAG = os.environ.get("SET_TAG", TAG)
PCS = int(os.environ.get("PCS", 0))
conds = pd.read_csv(CA / "z_conds.tsv", sep="\t").iloc[:, 0].tolist()
C = load.conditions().set_index("conditionID").loc[conds]
T = pd.read_csv(OUT / f"tables/specificity_per_kinase{TAG}.tsv", sep="\t", index_col=0)
kin = list(T.index)
spec = pd.read_parquet(CA / f"score_spec{SET_TAG}.parquet").loc[kin, conds]
conf = pd.read_parquet(CA / f"score_conf{SET_TAG}.parquet").loc[kin, conds]
full = pd.read_parquet(CA / f"score_full{SET_TAG}.parquet").loc[kin, conds]
b = T.beta_conf.fillna(0).values[:, None]
res = spec - b * conf.fillna(0).values
if PCS > 0:  # raw-space removal of global kinome modes: kinase loadings from study-centred full scores (03), applied to raw
    Fc = np.nan_to_num(pd.read_parquet(CA / f"full_sc{SET_TAG}.parquet").loc[kin, conds].values)
    Fc = Fc - Fc.mean(1, keepdims=True)
    U = np.linalg.svd(Fc, full_matrices=False)[0][:, :PCS]
    Fr = np.nan_to_num(full.values); Sr = U.T @ (Fr - Fr.mean(1, keepdims=True))  # PCS x conds
    R0 = res.values.copy()
    for i in range(len(kin)):
        y = spec.values[i]; m = ~np.isnan(y)
        if m.sum() < 30: continue
        X = np.c_[np.ones(len(conds)), np.nan_to_num(conf.values[i]), Sr.T]
        bb = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
        R0[i] = y - X @ bb + bb[0]
    res = pd.DataFrame(R0, index=kin, columns=conds)
sign = np.where(C.contrast_reversed.astype(bool).values, -1, 1)


def rstd(D):
    D = D * sign
    s = 1.4826 * (D.sub(D.median(1), axis=0)).abs().median(1)
    return D.div(s, axis=0)


S = {"resid": rstd(res), "conf": rstd(conf), "full": rstd(full)}
# context map
grp = {"MEK": ["MAP2K1", "MAP2K2"], "ERK": ["MAPK1", "MAPK3"], "RAF": ["BRAF", "ARAF", "RAF1"], "EGFR/ERBB": ["EGFR", "ERBB2"],
       "MET": ["MET"], "ALK/ROS1": ["ALK", "ROS1"], "MTOR": ["MTOR"], "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
       "AKT": ["AKT1", "AKT2", "AKT3"], "ATR/ATM/DNAPK": ["ATR", "ATM", "PRKDC"], "CDK": ["CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CDK12"],
       "p38": ["MAPK14", "MAPK11"], "JNK": ["MAPK8", "MAPK9", "MAPK10"], "PKC": ["PRKCA", "PRKCB", "PRKCG", "PRKCD", "PRKCE"],
       "GSK3": ["GSK3B", "GSK3A"], "ABL/SRC": ["ABL1", "SRC", "YES1", "LCK", "FYN"], "KIT/PDGFR": ["KIT", "PDGFRA", "PDGFRB"],
       "AURK/PLK": ["AURKA", "AURKB", "PLK1"], "JAK": ["JAK1", "JAK2", "JAK3"], "IKK/TBK1": ["IKBKB", "TBK1", "IKBKE"], "CK2": ["CSNK2A1", "CSNK2A2"]}
inv = {t: g for g, ts in grp.items() for t in ts}
ctx = {}
for c, r in C.iterrows():
    L = []
    if r.perturbation_class == "kinase_inhibitor" and isinstance(r.inhibitor_primary_targets, str) and r.perturbation_subclass in ("inhibitor", "inhibitor+stimulus"):
        L += sorted({"inh:" + inv[t] for t in r.inhibitor_primary_targets.split(";") if t in inv})
    elif r.perturbation_class != "kinase_inhibitor":
        if r.perturbation_subclass in ("antibody_CD20",) or r.perturbation_class == "time_course": pass
        else:
            L.append(f"{r.perturbation_class}:{r.perturbation_subclass}")
            if r.perturbation_class in ("growth_factor", "cytokine") and isinstance(r.stimulus, str):
                L.append(f"stim:{r.stimulus}")
    ctx[c] = L
cl = pd.Series(ctx).explode().dropna()
study = C.study_pmid
rows = []
for cx, cs in cl.groupby(cl):
    cs = cs.index.tolist()
    st = study.loc[cs]
    if st.nunique() < 2: continue
    for typ, D in S.items():
        sub = D[cs]
        sm = sub.T.groupby(st.values).mean().T  # kinase x study
        nst = sm.notna().sum(axis=1)
        mu = sm.mean(axis=1)
        for k in kin:
            v = sm.loc[k].dropna().values
            if len(v) < 2: continue
            p = ttest_1samp(v, 0).pvalue if len(v) >= 3 else np.nan
            rows.append(dict(kinase=k, context=cx, score=typ, n_cond=int(sub.loc[k].notna().sum()), n_study=len(v),
                             effect=v.mean(), frac_same_sign=np.mean(np.sign(v) == np.sign(v.mean())), p=p))
R = pd.DataFrame(rows)
R["fdr"] = np.nan
for typ in S:
    m = (R.score == typ) & R.p.notna()
    R.loc[m, "fdr"] = multipletests(R.loc[m, "p"], method="fdr_bh")[1]
R["dark"] = T.dark.reindex(R.kinase).values
R.to_csv(OUT / f"tables/contexts_long{TAG}.tsv", sep="\t", index=False)
W = R.pivot_table(index=["kinase", "context"], columns="score", values=["effect", "p", "fdr"]).reset_index()
W.columns = ["_".join([c for c in col if c]) for col in W.columns]
W = W.merge(R[R.score == "resid"][["kinase", "context", "n_study", "n_cond", "frac_same_sign"]], on=["kinase", "context"])
W["dark"] = T.dark.reindex(W.kinase).values
W.to_csv(OUT / f"tables/contexts_wide{TAG}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
hits = W[(W.fdr_resid <= 0.1) & W.dark].sort_values("p_resid")
print(len(hits), "dark kinase-context hits FDR<=0.1 (resid)")
print(hits.head(40).round(4).to_string())
print(R[R.score == "resid"].groupby("dark").apply(lambda d: (d.fdr <= 0.1).sum()))
