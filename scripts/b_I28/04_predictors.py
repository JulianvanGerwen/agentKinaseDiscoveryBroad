"""I28 step 4: what predicts fast turnover?
(a) per-feature logistic (fast vs stable), family fixed effects + coverage covariates, cluster-robust SE by protein;
(b) multivariable L2 logistic / ridge; (c) leave-one-study-out validation with sites ALSO held out (pure feature generalisation).
Out: tables/predictors_univariate.tsv, tables/predictors_multivariable.tsv, tables/predictors_loso.tsv
"""
import os, sys, warnings
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, LogisticRegression
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(__file__)); from lib28 import *
warnings.filterwarnings("ignore")
rng = np.random.default_rng(4)
A = pd.read_csv(f"{R28}/tables/turnover_atlas.tsv", sep="\t")
F = pd.read_parquet(f"{C28}/features.parquet").drop(columns=["EV_max_frequency", "isHotspot", "isInterface"])
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]
A = A[A.group.isin(MAIN) & (A.n_units >= 2)].drop(columns=["psp_regulatory"]).copy()
A["prot"] = A.site.str.split("_").str[0]
A["motif_pct_f"] = A.motif_pct; A["is_psp"] = (A.source != "pred").astype(float); A["log_psp_LT"] = np.log10(1 + A.psp_LT.fillna(0))
A["n_units_f"] = np.log2(A.n_units)
D = A.join(F, on="site")
FEATS = [c for c in F.columns if c not in ("log_n_conditions", "log_n_studies", "quant_top1")] + ["motif_pct_f", "is_psp", "log_psp_LT"]
COV = ["log_n_conditions", "n_units_f"]
grp = pd.get_dummies(D.group, drop_first=True, dtype=float)

# ---------- (a) univariate (each feature + group FE + coverage)
rows = []
fs = D[D.cls.isin(["fast", "stable"])]
for f in FEATS + ["quant_top1", "log_n_conditions", "lability_other"]:
    d = fs.dropna(subset=[f] + COV)
    if d[f].std() == 0 or len(d) < 50:
        continue
    x = (d[f] - d[f].mean()) / d[f].std()
    X = pd.concat([x.rename(f), grp.loc[d.index]] + ([d[COV]] if f not in COV else [d[[c for c in COV if c != f]]]), axis=1)
    X = sm.add_constant(X)
    y = (d.cls == "fast").astype(float)
    try:
        m = sm.Logit(y, X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": d.prot.astype("category").cat.codes})
        b, se, p = m.params[f], m.bse[f], m.pvalues[f]
    except Exception as e:
        continue
    # continuous: within-group Spearman of feature vs -R (all classes)
    dd = D.dropna(subset=[f])
    rr = [spearmanr(g[f], -g.R)[0] for _, g in dd.groupby("group") if len(g) > 20 and g[f].std() > 0]
    rows.append(dict(feature=f, n=len(d), OR_per_SD=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=p,
                     mean_fast=d[f][y == 1].mean(), mean_stable=d[f][y == 0].mean(), spearman_negR_within_group=np.nanmean(rr)))
U = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
U["q"] = multipletests(U.p, method="fdr_bh")[1]
U = U.sort_values("p"); U.to_csv(f"{R28}/tables/predictors_univariate.tsv", sep="\t", index=False)
print(U.round(3).to_string())


# ---------- (b) multivariable
def design(d, feats, med=None, sd=None):
    X = d[feats].astype(float)
    miss = X.isna().astype(float).add_suffix("_na"); miss = miss.loc[:, [c for c in miss.columns if c.replace("_na", "") in ("functional_score", "disopred_score", "w0_mya", "adj_ptms_w21", "paxdb_abundance_log10", "motif_pct_f")]]
    if med is None:
        med = X.median(); sd = X.std().replace(0, 1)
    X = ((X.fillna(med) - med) / sd)
    g = pd.get_dummies(d.group, dtype=float).reindex(columns=MAIN, fill_value=0)
    return pd.concat([X, miss, g], axis=1), med, sd


INTR = [f for f in FEATS]            # intrinsic + annotation features
fsx = D[D.cls.isin(["fast", "stable"])]
X, med, sd = design(fsx, INTR + COV)
y = (fsx.cls == "fast").astype(int)
lr = LogisticRegression(C=0.3, max_iter=5000).fit(X, y)
coef = pd.Series(lr.coef_[0], index=X.columns)
# protein-cluster bootstrap CIs
prots = fsx.prot.unique(); bsc = []
for b in range(200):
    pb = rng.choice(prots, len(prots)); idx = np.concatenate([np.where(fsx.prot.values == p)[0] for p in pb])
    if y.iloc[idx].nunique() < 2: continue
    bsc.append(LogisticRegression(C=0.3, max_iter=3000).fit(X.iloc[idx], y.iloc[idx]).coef_[0])
bsc = np.array(bsc)
MV = pd.DataFrame({"coef_logOR_per_SD": coef, "lo": np.percentile(bsc, 2.5, 0), "hi": np.percentile(bsc, 97.5, 0)})
MV["sig"] = (MV.lo > 0) | (MV.hi < 0)
MV.sort_values("coef_logOR_per_SD").to_csv(f"{R28}/tables/predictors_multivariable.tsv", sep="\t")
print(MV.sort_values("coef_logOR_per_SD").round(3).to_string())

# ---------- (c) LOSO with held-out sites
L = pd.read_parquet(f"{C28}/L_units.parquet"); L = L[(L.set == "substrate") & L.group.isin(MAIN)]
Fj = F
res = []
for S_ in sorted(L.study.unique()):
    te = L[L.study == S_].groupby(["group", "site"]).z.mean().rename("zS").reset_index()
    tr = classify_units(L[L.study != S_])
    tr = tr[tr.n_units >= 1]
    tr_other = tr.set_index(["group", "site"]).R
    te["R_other"] = tr_other.reindex(pd.MultiIndex.from_frame(te[["group", "site"]])).values
    trn = tr[~tr.site.isin(set(te.site))].merge(A[["group", "site", "motif_pct_f", "is_psp", "log_psp_LT", "lability_other"]], on=["group", "site"], how="left")
    trn["n_units_f"] = np.log2(trn.n_units)
    trn = trn.join(Fj, on="site")
    if len(te) < 30 or trn.group.nunique() < 1:
        continue
    teD = te.merge(A[["group", "site", "motif_pct_f", "is_psp", "log_psp_LT", "lability_other"]], on=["group", "site"], how="left").join(Fj, on="site")
    teD = teD[teD.motif_pct_f.notna() | teD.is_psp.notna()]
    teD["n_units_f"] = 0.0
    out = {}
    for name, fl in [("features", INTR + ["log_n_conditions"]), ("motif_only", ["motif_pct_f"]), ("features+lability", INTR + ["log_n_conditions", "lability_other"])]:
        Xtr, m_, s_ = design(trn, fl); Xte, _, _ = design(teD, fl, m_, s_)
        Xte = Xte.reindex(columns=Xtr.columns, fill_value=0)
        w = np.clip(trn.n_units.values, 1, 4)
        mdl = RidgeCV(alphas=np.logspace(0, 4, 9)).fit(Xtr, -trn.R.clip(-6, 4), sample_weight=w)
        out[name] = mdl.predict(Xte)
    for name, p in out.items():
        teD[f"pred_{name}"] = p
    teD["study"] = S_
    res.append(teD[["group", "site", "study", "zS", "R_other"] + [f"pred_{n}" for n in out]])
    print(S_, len(teD), flush=True)
P = pd.concat(res)
P.to_parquet(f"{C28}/loso_predictions.parquet")


def within(P, col):
    out = []
    for (s, g), d in P.groupby(["study", "group"]):
        d = d.dropna(subset=[col, "zS"])
        if len(d) >= 20 and d[col].std() > 0:
            out.append((s, g, len(d), spearmanr(d[col], -d.zS)[0], roc_auc_score((d.zS <= -1).astype(int), d[col]) if 0 < (d.zS <= -1).sum() < len(d) else np.nan))
    return pd.DataFrame(out, columns=["study", "group", "n", "rho", "auroc_fall"])


rows = []
for col in ["pred_features", "pred_motif_only", "pred_features+lability", "R_other"]:
    x = P.copy()
    if col == "R_other":
        x[col] = -x[col]
    w = within(x, col); w["predictor"] = col; rows.append(w)
W = pd.concat(rows)
W.to_csv(f"{R28}/tables/predictors_loso.tsv", sep="\t", index=False)
summ = W.groupby("predictor").apply(lambda d: pd.Series(dict(n_strata=len(d), n_sites=d.n.sum(), rho_wmean=np.average(d.rho, weights=d.n),
                                                              auroc_wmean=np.average(d.auroc_fall.fillna(0.5), weights=d.n))))
print(W.round(3).to_string()); print(summ.round(3).to_string())
summ.to_csv(f"{R28}/tables/predictors_loso_summary.tsv", sep="\t")
