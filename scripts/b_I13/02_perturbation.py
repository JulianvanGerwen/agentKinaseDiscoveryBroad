"""Task 1: which sequence features predict increases under phosphatase inhibitors (compendium)."""
import sys, json
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from sklearn.inspection import permutation_importance
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
R = "/workspace/results/b_I13"
F = pd.read_parquet(f"{R}/cache/site_features.parquet").set_index("site")
S = load.sites().set_index("site"); C = load.conditions()
X = load.sites_matrix()
# generic site dynamism: robust SD of log2FC across all conditions (>=20 measurements)
n = X.notna().sum(1)
med = X.median(1); mad = (X.sub(med, axis=0)).abs().median(1)
dyn = np.log10(1.4826 * mad + 0.05).where(n >= 20)
F["dyn"] = dyn.reindex(F.index); F["log_nmeas"] = np.log10(n.reindex(F.index).fillna(0) + 1)
F.loc[F.dyn.isna(), "dyn"] = F.dyn.median()

def rz(v):
    v = v.dropna(); return (v - v.median()) / (1.4826 * (v - v.median()).abs().median())
PHOS = ["isT", "isY", "P1", "Q1", "R_m3", "R_m2", "basic_m7m1", "basic_p1p7", "bipartite_basic", "acid_m3p3", "acid_p1p3",
        "pro_m7p7", "hydrophob_p1", "ST_m7p7", "logdist_LxxIxE", "logdist_RVxF", "logdist_PxIxIT", "logdist_LxVP", "prot_len_log"]
NUIS = ["dyn", "log_nmeas"]
# kinase-motif baseline: multi-hot of top-3 predicted kinases (kinaseSubMain) for frequent kinases
pk = S.pred_top3_kinases.fillna("").str.split(";")
from collections import Counter
cnt = Counter(k for l in pk for k in l if k)
KIN = [k for k, c in cnt.items() if c >= 500]
KM = pd.DataFrame({f"kin_{k}": pk.apply(lambda l, k=k: k in l).astype(int) for k in KIN}, index=S.index)
F = F.join(KM, how="left"); F[KM.columns] = F[KM.columns].fillna(0)
KINF = list(KM.columns)

SETS = {"PPP_inhib": {"OA": "26060313_Ochoa_1306_272_log2", "CA": "26060313_Ochoa_1307_272_log2"},
        "PTP_pervanadate": {"PV_AML_23532336": "23532336_qPTM_6_log2", "PV_GMP_23532336": "23532336_qPTM_7_log2", "PV_HeLa_28674419": "28674419_qPTM_1_log2"}}
rows, coefs = [], []
def fit_eval(tr, te, feats, kind="lr"):
    if kind == "lr":
        sc = StandardScaler().fit(tr[feats]); m = LogisticRegression(C=0.5, max_iter=2000, class_weight="balanced").fit(sc.transform(tr[feats]), tr.y)
        p = m.predict_proba(sc.transform(te[feats]))[:, 1]
    else:
        m = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=200, class_weight="balanced", random_state=0).fit(tr[feats], tr.y)
        p = m.predict_proba(te[feats])[:, 1]
    return roc_auc_score(te.y, p), m
for setname, conds in SETS.items():
    D = []
    for lab, cid in conds.items():
        z = rz(X[cid]); d = F.reindex(z.index).dropna(subset=["isT"]).copy(); d["z"] = z.reindex(d.index); d["y"] = (d.z >= 2).astype(int)
        d["cond"] = lab; d["study"] = cid.split("_")[0]; D.append(d)
    D = pd.concat(D)
    FS = {"nuisance": NUIS, "kinase_motif": NUIS + KINF, "phosphatase_feats": NUIS + PHOS, "kinase+phosphatase": NUIS + KINF + PHOS}
    # held-out: leave-one-condition-out (PPP: same study, 2 drugs); PTP: leave-one-study-out
    grp = "cond" if setname == "PPP_inhib" else "study"
    for g in D[grp].unique():
        tr, te = D[D[grp] != g], D[D[grp] == g]
        for fsn, fs in FS.items():
            for kind in ["lr", "gbm"]:
                auc, _ = fit_eval(tr, te, fs, kind)
                rows.append(dict(set=setname, heldout=g, features=fsn, model=kind, auc=auc, n_test=len(te), n_pos=int(te.y.sum())))
    # univariate per-condition effects (logistic, nuisance-adjusted), standardized
    import statsmodels.api as sm
    for lab in D.cond.unique():
        d = D[D.cond == lab]
        for f in PHOS:
            if d[f].std() == 0: continue
            Xm = sm.add_constant(pd.DataFrame({"f": (d[f] - d[f].mean()) / d[f].std(), "dyn": d.dyn, "nm": d.log_nmeas}))
            try:
                r = sm.Logit(d.y, Xm).fit(disp=0); coefs.append(dict(set=setname, cond=lab, feature=f, beta=r.params["f"], se=r.bse["f"], p=r.pvalues["f"]))
            except Exception: pass
    # permutation importance of final GBM on pooled (with GroupKFold by protein for honest importance)
    m = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=200, class_weight="balanced", random_state=0)
    feats = NUIS + KINF + PHOS
    gkf = GroupKFold(5); imps = []
    for tri, tei in gkf.split(D, D.y, D.uniprot):
        m.fit(D.iloc[tri][feats], D.iloc[tri].y)
        pi = permutation_importance(m, D.iloc[tei][feats], D.iloc[tei].y, scoring="roc_auc", n_repeats=3, random_state=0, n_jobs=1)
        imps.append(pi.importances_mean)
    pd.DataFrame({"feature": feats, "perm_importance_auc": np.mean(imps, 0)}).sort_values("perm_importance_auc", ascending=False).to_csv(f"{R}/tables/t1_{setname}_gbm_perm_importance.tsv", sep="\t", index=False)
pd.DataFrame(rows).to_csv(f"{R}/tables/t1_heldout_auc.tsv", sep="\t", index=False)
pd.DataFrame(coefs).to_csv(f"{R}/tables/t1_univariate_logit.tsv", sep="\t", index=False)

# specificity control: per-condition pT-vs-pS and bipartite-basic effect across ALL compendium conditions (rank of OA/CA)
Fx = F.reindex(X.index)
ctrl = []
for cid in C.conditionID[(C.qc_flag == "ok")]:
    z = rz(X[cid]);
    if len(z) < 2000: continue
    f = Fx.reindex(z.index); st = f.residue.isin(["S", "T"])
    zz, ff = z[st], f[st]
    ctrl.append(dict(cond=cid, cls=C.set_index("conditionID").perturbation_class[cid],
                     dT=zz[ff.isT == 1].mean() - zz[ff.isT == 0].mean(),
                     dP1=zz[ff.P1 == 1].mean() - zz[ff.P1 == 0].mean(),
                     dBip=zz[ff.bipartite_basic == 1].mean() - zz[ff.bipartite_basic == 0].mean(),
                     dB55=np.corrcoef(zz, ff.B55_lit)[0, 1],
                     dY=z[f.isY == 1].mean() - z[f.isY == 0].mean()))
ctrl = pd.DataFrame(ctrl)
for c in ["dT", "dP1", "dBip", "dB55", "dY"]: ctrl[c + "_pct"] = ctrl[c].rank(pct=True)
ctrl.to_csv(f"{R}/tables/t1_specificity_all_conditions.tsv", sep="\t", index=False)
allc = sum([list(v.values()) for v in SETS.values()], [])
print(ctrl[ctrl.cond.isin(allc)].round(3).to_string())
print(pd.DataFrame(rows).pivot_table(index=["set", "features"], columns="model", values="auc", aggfunc="mean").round(3))
