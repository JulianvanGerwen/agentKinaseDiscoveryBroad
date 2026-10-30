"""Kinase 'reach': fraction of the measured phosphoproteome responding per unit target occupancy.
Occupancy at the site EC50: theta = EC50/(EC50+Kd_target) (Klaeger Kd_app, simple 1:1 binding).
reach_theta50 = frac measured sites with EC50 <= Kd (respond by 50% occupancy); reach_theta90 = EC50 <= 9 Kd;
reach_any = frac significant at any dose. Related to target centrality (STRING>=700 degree, OmniPath out-degree,
PSP substrate count) and DepMap mean Chronos (pan-line; no model map available on disk).
Output: results/b_I17/reach_by_experiment.tsv, reach_correlations.tsv
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"
F = pd.read_parquet(f"{RES}/curves_all.parquet")
X = pd.read_csv(f"{RES}/shift_by_experiment.tsv", sep="\t")
ST = pd.read_parquet(f"{OUT}/string_edges.parquet")
deg = pd.concat([ST.a, ST.b]).value_counts() / 2
om = pd.read_parquet(f"{OUT}/omnipath_protein_edges.parquet")
outdeg = om.groupby("a").b.nunique()
KA = pd.read_csv("/workspace/results/infra_compendium/kinase_annotation.tsv", sep="\t")
dm = pd.read_csv(f"{OUT}/depmap_mean.tsv", sep="\t", index_col=0)
rows = []
for _, x in X.iterrows():
    f = F[F.exp_id == x.exp_id]
    kd = x.kd_target_nM
    ec = 10 ** f.logEC50
    sig = f.sig
    rows.append(dict(exp_id=x.exp_id, agent=x.agent, cell_line=x.cell_line, target=x.target, kd_nM=kd, n_measured=len(f),
                     reach_any=sig.mean(), reach_theta50=(sig & (ec <= kd)).mean(), reach_theta90=(sig & (ec <= 9 * kd)).mean(),
                     reach_down_theta90=(sig & (f.delta < 0) & (ec <= 9 * kd)).mean(),
                     string_degree=deg.get(x.target, np.nan), omnipath_outdeg=outdeg.get(x.target, np.nan),
                     mean_chronos=dm.mean_chronos.get(x.target, np.nan), frac_dep=dm.frac_dep.get(x.target, np.nan)))
R = pd.DataFrame(rows); R.to_csv(f"{RES}/reach_by_experiment.tsv", sep="\t", index=False)
print(R.round(4).to_string())
# one row per target (median over experiments) to avoid pseudo-replication
T = R.groupby("target").median(numeric_only=True)
cr = []
for y in ["reach_any", "reach_theta90"]:
    for x in ["string_degree", "omnipath_outdeg", "mean_chronos"]:
        r, p = spearmanr(T[x], T[y], nan_policy="omit")
        cr.append(dict(reach=y, covariate=x, rho=r, p=p, n_targets=int(T[[x, y]].dropna().shape[0])))
C = pd.DataFrame(cr); C.to_csv(f"{RES}/reach_correlations.tsv", sep="\t", index=False); print(T.round(4).to_string()); print(C.round(3).to_string())
