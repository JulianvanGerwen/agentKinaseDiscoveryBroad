"""Baseline (DepMap target expression) vs reverse-layer LR: added value; LODO CV. Out: tables/T8_baseline_compare.tsv"""
import os
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]: os.environ[v]="4"
import sys; sys.path.insert(0, "/workspace/scripts/b_I23")
import numpy as np, pandas as pd, re
from importlib import import_module
P = import_module("04_primary")
T = "/workspace/results/b_I23/tables/"; E = "/workspace/cache/b_I23/dl_expr/"
D = P.D.copy()
M = pd.read_csv(E + "Model.csv")
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())
M["sn"] = M.StrippedCellLineName.map(norm)
alias = {"NTERA2": "NTERA2CLD1"}
ids = {c: M.loc[M.sn == alias.get(c, c), "ModelID"].tolist() for c in D.cell.unique()}
print(ids)
X = pd.read_csv(E + "expr.csv", index_col=0)
X.columns = [c.split(" (")[0] for c in X.columns]
Xz = (X - X.mean()) / X.std()
Dt = pd.read_csv("/workspace/results/infra_compendium/drug_targets.tsv", sep="\t").set_index("agent")
def texpr(r):
    tg = [t for t in str(Dt.primary_targets_in_universe.get(r.drug, "")).split(";") if t in Xz.columns]
    mid = [i for i in ids[r.cell] if i in Xz.index]
    return Xz.loc[mid[0], tg].mean() if tg and mid else np.nan
D["target_expr"] = D.apply(texpr, axis=1)
print("target_expr non-missing", D.target_expr.notna().sum())
R = []
R.append(dict(test="baseline_target_expr", **P.run(D, "auc", ["target_expr"], nperm=2000, nboot=500)))
R.append(dict(test="LR_adj_target_expr", **P.run(D, "auc", ["LR", "target_expr"], nperm=2000, nboot=500)))
R.append(dict(test="LR_adj_expr_targetz", **P.run(D, "auc", ["LR", "target_expr", "target_z"], nperm=2000, nboot=500)))
# LODO CV on FE residuals
def cv(df, xs):
    df = df.dropna(subset=["auc"] + xs).reset_index(drop=True)
    Y = P.resid(df, "auc"); Xm = np.column_stack([P.resid(df, x) for x in xs])
    pred = np.zeros(len(df))
    for d in df.drug.unique():
        te = (df.drug == d).values
        b, *_ = np.linalg.lstsq(Xm[~te], Y[~te], rcond=None); pred[te] = Xm[te] @ b
    return 1 - ((Y - pred) ** 2).sum() / (Y ** 2).sum(), np.corrcoef(Y, pred)[0, 1], len(df)
sub = D.dropna(subset=["target_expr", "target_z"])
for xs in [["target_expr"], ["target_z"], ["LR"], ["LR", "target_expr"], ["LR", "target_expr", "target_z"], ["logfu"], ["logfd"]]:
    r2, rr, n = cv(sub, xs); R.append(dict(test="LODO_CV_" + "+".join(xs), n=n, cv_R2=r2, cv_r=rr))
R = pd.DataFrame(R); R.to_csv(T + "T8_baseline_compare.tsv", sep="\t", index=False)
D.to_csv(T + "T3_hijazi_pairs.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(R.round(4).to_string())
