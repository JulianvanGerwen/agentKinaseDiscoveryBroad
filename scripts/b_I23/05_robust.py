"""Robustness of P1 + S3 kinase groups. Out: tables/T5_robustness.tsv, T6_kinase_groups.tsv, T7_per_drug.tsv"""
import os
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]: os.environ[v]="4"
import sys; sys.path.insert(0, "/workspace/scripts/b_I23")
import numpy as np, pandas as pd
from importlib import import_module
P = import_module("04_primary")  # reuses D and run()
from scipy import stats
T = "/workspace/results/b_I23/tables/"
D = P.D.copy()
F = pd.read_csv(T + "T1_features.tsv", sep="\t")
D["total"] = np.log(D.f_up + D.f_down); D["logmeas"] = np.log(D.n_meas)
R = []
R.append(dict(test="adj_total_magnitude", **P.run(D, "auc", ["LR", "total"], nperm=2000, nboot=500)))
R.append(dict(test="adj_n_measured", **P.run(D, "auc", ["LR", "logmeas"], nperm=2000, nboot=500)))
R.append(dict(test="GDSC2_only", **P.run(D[D.source == "GDSC2"], "auc", ["LR"], nperm=2000, nboot=500)))
R.append(dict(test="GDSC1_only", **P.run(D[D.source == "GDSC1"], "auc", ["LR"], nperm=2000, nboot=500)))
for c in ["HL60", "MCF7", "NTERA2"]:
    d = D[D.cell != c]; d = d[d.groupby("drug").cell.transform("size") == 2]
    R.append(dict(test=f"drop_{c}", **P.run(d, "auc", ["LR"], nperm=2000, nboot=500)))
# rank outcome within drug
D["auc_rank"] = D.groupby("drug").auc.rank()
R.append(dict(test="rank_outcome", **P.run(D, "auc_rank", ["LR"], nperm=2000, nboot=500)))
# shuffled exposure negative control (one draw of across-everything shuffle, 200 reps summarised)
rng = np.random.default_rng(5); sb = []
for i in range(200):
    d = D.copy(); d["LRs"] = rng.permutation(d.LR.values)
    b, Y, X = P.fe_beta(d, "auc", ["LRs"]); sb.append(b[0] * X[:, 0].std() / Y.std())
R.append(dict(test="shuffled_LR_null_mean_sd", std_beta=np.mean(sb), ci_lo=np.percentile(sb, 2.5), ci_hi=np.percentile(sb, 97.5)))
# LODO
lodo = []
for d0 in D.drug.unique():
    d = D[D.drug != d0]; b, Y, X = P.fe_beta(d, "auc", ["LR"]); lodo.append(b[0] * X[:, 0].std() / Y.std())
R.append(dict(test="LODO_range", std_beta=np.median(lodo), ci_lo=min(lodo), ci_hi=max(lodo)))
R = pd.DataFrame(R); R.to_csv(T + "T5_robustness.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(R.round(4).to_string())
# per-drug
pd_rows = []
for d0, g in D.groupby("drug"):
    pd_rows.append(dict(drug=d0, n=len(g), source=g.source.iloc[0], auc=";".join(f"{c}:{a:.2f}" for c, a in zip(g.cell, g.auc)),
                        LR=";".join(f"{c}:{a:.2f}" for c, a in zip(g.cell, g.LR)), spearman=stats.spearmanr(g.LR, g.auc)[0] if len(g) > 2 else np.sign(np.corrcoef(g.LR, g.auc)[0, 1]),
                        auc_range=g.auc.max() - g.auc.min()))
PD = pd.DataFrame(pd_rows).sort_values("auc_range", ascending=False); PD.to_csv(T + "T7_per_drug.tsv", sep="\t", index=False)
print(PD.to_string()); print("sign counts", np.sign(PD.spearman).value_counts().to_dict())
# S3 kinase groups
G = ["STE20_MAP4K", "NEK", "IRAK4", "MAP3K", "MARK_SIK", "AKT", "SFK"]
K = []
for g in G:
    for pre in ["z_", "motif_"]:
        col = pre + g
        if col not in D or D[col].notna().sum() < 40: continue
        r = P.run(D, "auc", [col], nperm=2000, nboot=500)
        r2 = P.run(D, "auc", [col, "LR"], nperm=1000, nboot=200)
        lo = []
        dd = D.dropna(subset=[col])
        for d0 in dd.drug.unique():
            q = dd[dd.drug != d0]; b, Y, X = P.fe_beta(q, "auc", [col]); lo.append(np.sign(b[0]))
        K.append(dict(group=g, feature=pre.strip("_"), n=r["n"], std_beta=r["std_beta"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], p=r["p_perm"],
                      std_beta_adjLR=r2["std_beta"], p_adjLR=r2["p_perm"], lodo_sign_frac=np.mean(np.array(lo) == np.sign(r["beta"]))))
K = pd.DataFrame(K)
from statsmodels.stats.multitest import multipletests
for f in K.feature.unique():
    m = K.feature == f; K.loc[m, "q_BH"] = multipletests(K.loc[m, "p"], method="fdr_bh")[1]
K.to_csv(T + "T6_kinase_groups.tsv", sep="\t", index=False); print(K.round(3).to_string())
