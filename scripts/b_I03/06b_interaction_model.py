"""I03 step 6b: does the new evidence add to the sequence/PPI prior at the TOP of the list?
Gradient boosting (monotone constraints, shallow) trained on DEV proteins (PSP labels), tested on held-out TEST proteins
(PSP) and on curated external non-PSP edges. Prior = kinaseSub DEV/no-phospho pct + Cantley motif + PPI;
FULL = prior + perturbation + dose + tumour. Metrics: macro AUROC, precision at top-N, top-1% enrichment.
Also writes stratified PSP-rate table (kinaseSub tier x replicated drop x tumour).
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
OUT = "/workspace/cache/b_I03"; RES = "/workspace/results/b_I03/tables"
X = pd.read_parquet(f"{OUT}/scored.parquet")
rng = np.random.default_rng(0)
PRIOR = ["f_ksdev", "f_motif", "f_ppi"]
NEWF = ["f_rep", "f_stouf", "f_ndrop", "f_spec", "f_up", "f_dose", "f_dosehit", "f_tum"]
mono = {"f_ksdev": 1, "f_motif": 1, "f_ppi": 1, "f_rep": 1, "f_stouf": 1, "f_ndrop": 1, "f_spec": 0, "f_up": 0, "f_dose": 1, "f_dosehit": 1, "f_tum": 1}
def fit(feats, tr, y, seed=0):
    m = HistGradientBoostingClassifier(max_depth=3, max_iter=300, learning_rate=0.05, l2_regularization=1.0, min_samples_leaf=200,
                                       monotonic_cst=[mono[f] for f in feats], class_weight="balanced", random_state=seed)
    m.fit(X.loc[tr, feats].values, y[tr]); return m
def aucs(s, y, kin):
    out = []
    for k, idx in pd.Series(np.arange(len(y))).groupby(kin.values).groups.items():
        if y[idx].sum() >= 3: out.append(roc_auc_score(y[idx], s[idx]))
    return np.mean(out), len(out)
def met(s, y, kin):
    o = np.argsort(-s); ma, nk = aucs(s, y, kin)
    return dict(macro_auroc=ma, n_kin=nk, prec_top100=y[o[:100]].mean(), prec_top300=y[o[:300]].mean(),
                enrich_top1pct=y[o[:max(1, len(y) // 100)]].mean() / y.mean(), n_pos=int(y.sum()), n=len(y))
y = X.y_psp.values
tr = (X.test_protein == 0).values
rows = []
for nm, feats in [("PRIOR (kinaseSub-DEV + motif + PPI)", PRIOR), ("FULL (PRIOR + perturbation + dose + tumour)", PRIOR + NEWF),
                  ("NEW only", NEWF)]:
    ss = []
    for seed in range(3):
        m = fit(feats, tr, y, seed); ss.append(m.predict_proba(X[feats].values)[:, 1])
    X["g_" + nm.split()[0]] = np.mean(ss, 0)
te = X.test_protein == 1
T = X[te]; E = X[(X.y_psp == 0)]
for nm in ["PRIOR", "FULL", "NEW"]:
    rows.append(dict(eval="PSP TEST proteins", model=nm, **met(T["g_" + nm].values, T.y_psp.values, T.kinase)))
    rows.append(dict(eval="external curated non-PSP (all proteins)", model=nm, **met(E["g_" + nm].values, E.y_ext.values, E.kinase)))
    Et = E[E.test_protein == 1]
    rows.append(dict(eval="external curated non-PSP (TEST proteins)", model=nm, **met(Et["g_" + nm].values, Et.y_ext.values, Et.kinase)))
# permutation control: shuffle the new-evidence rows across sites within kinase (breaks site link, keeps distributions)
Xp = X.copy()
for f in NEWF:
    Xp[f] = Xp.groupby("kinase")[f].transform(lambda v: v.sample(frac=1, random_state=1).values)
Xs = X; X = Xp
m = fit(PRIOR + NEWF, tr, y, 0); X = Xs
X["g_FULLperm"] = m.predict_proba(Xp[PRIOR + NEWF].values)[:, 1]
T = X[te]; E = X[(X.y_psp == 0)]
rows.append(dict(eval="PSP TEST proteins", model="FULL, new evidence shuffled within kinase", **met(T.g_FULLperm.values, T.y_psp.values, T.kinase)))
rows.append(dict(eval="external curated non-PSP (all proteins)", model="FULL, new evidence shuffled within kinase", **met(E.g_FULLperm.values, E.y_ext.values, E.kinase)))
R = pd.DataFrame(rows); R.to_csv(f"{RES}/val_gbm_prior_vs_full.tsv", sep="\t", index=False); print(R.round(3).to_string())
# bootstrap over kinases for FULL - PRIOR
def bd(D, ycol):
    a = []; 
    for k, d in D.groupby("kinase"):
        if d[ycol].sum() >= 3: a.append(roc_auc_score(d[ycol], d.g_FULL) - roc_auc_score(d[ycol], d.g_PRIOR))
    a = np.array(a); bs = [a[rng.integers(0, len(a), len(a))].mean() for _ in range(1000)]
    return a.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5)
print("FULL-PRIOR TEST", bd(T, "y_psp"), "external", bd(E, "y_ext"))
# precision-at-N bootstrap (sites resampled) for TEST
def pdiff(D, ycol, N=200, B=500):
    out = []
    for _ in range(B):
        d = D.sample(frac=1, replace=True, random_state=int(rng.integers(1e9)))
        out.append(d.nlargest(N, "g_FULL")[ycol].mean() - d.nlargest(N, "g_PRIOR")[ycol].mean())
    return np.mean(out), np.percentile(out, 2.5), np.percentile(out, 97.5)
pt = pdiff(T, "y_psp"); pe = pdiff(E, "y_ext")
print("prec@200 FULL-PRIOR TEST", pt, "external", pe)
pd.DataFrame([dict(eval="PSP TEST", stat="macro AUROC FULL-PRIOR", **dict(zip(["mean", "lo", "hi"], bd(T, "y_psp")))),
              dict(eval="external", stat="macro AUROC FULL-PRIOR", **dict(zip(["mean", "lo", "hi"], bd(E, "y_ext")))),
              dict(eval="PSP TEST", stat="precision@200 FULL-PRIOR", **dict(zip(["mean", "lo", "hi"], pt))),
              dict(eval="external", stat="precision@200 FULL-PRIOR", **dict(zip(["mean", "lo", "hi"], pe)))]).to_csv(f"{RES}/val_gbm_bootstrap.tsv", sep="\t", index=False)
# stratified table
X["ks_tier"] = pd.cut(X.ks_dev_pct.fillna(0), [-1, 95, 99, 99.7, 100.1], labels=["<95", "95-99", "99-99.7", ">99.7"])
X["pert"] = np.select([X.n_st_drop >= 2, X.n_st_drop == 1], ["drop in >=2 studies", "drop in 1 study"], "no drop")
X["tumour"] = np.select([X.tum_z > 3, X.tum_z.notna()], ["co-varies (z>3)", "no"], "not measured")
ST = X.groupby(["ks_tier", "pert", "tumour"], observed=True).agg(n=("y_psp", "size"), psp_rate=("y_psp", "mean"),
                                                                  ext_rate=("y_ext", "mean")).reset_index()
ST.to_csv(f"{RES}/val_stratified_psp_rate.tsv", sep="\t", index=False)
ST2 = X[X.test_protein == 1].groupby(["ks_tier", "pert"], observed=True).agg(n=("y_psp", "size"), psp_rate=("y_psp", "mean")).reset_index()
ST2.to_csv(f"{RES}/val_stratified_psp_rate_TEST.tsv", sep="\t", index=False); print(ST2.to_string())
X[["kinase", "site", "g_PRIOR", "g_FULL", "g_NEW", "g_FULLperm"]].to_parquet(f"{OUT}/gbm_scores.parquet")
