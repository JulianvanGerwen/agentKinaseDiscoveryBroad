"""I03 step 6: validation of single evidence lines and integrated scores.
Universe: panel kinase x site pairs measured in >=2 independent studies of K-loss (n_st>=2), residue-class matched.
Labels: PSP (label-free features -> fully held out) and curated external non-PSP edges (OmniPath non-PSP sources + iPTMnet).
Models (logistic regression, standardised features, substrate-protein grouped 5-fold CV):
  NEW        perturbation + dose + tumour (no sequence, no PPI, no kinaseSub)
  NEW+motif  + Cantley percentile
  integrated NEW+motif logit + kinaseSub pct + PPI, stacked on EXTERNAL edges (grouped OOF) -- kinaseSub is in-sample on PSP.
Outputs results/b_I03/tables/val_*.tsv, cache/b_I03/scored.parquet
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
OUT = "/workspace/cache/b_I03"; RES = "/workspace/results/b_I03/tables"
os.makedirs(RES, exist_ok=True)
rng = np.random.default_rng(0)
X = pd.read_parquet(f"{OUT}/pairs.parquet")
X = X[X.n_st >= 2].reset_index(drop=True)
X["uniprot"] = X.site.str.split("_").str[0]
# features
X["f_rep"] = -X.min2_z.clip(-8, 8)                  # replicated drop (2nd strongest study decrease)
X["f_stouf"] = -X.stouffer.clip(-10, 10)
X["f_ndrop"] = X.n_st_drop / X.n_st
X["f_spec"] = -(X.mean_z - X.gen_z.fillna(0)).clip(-8, 8)    # drop beyond the site's generic inhibitor response
X["f_up"] = X.n_st_up / X.n_st                      # reverse responses (penalty)
X["f_dose"] = X.dose_r.fillna(0); X["m_dose"] = X.dose_r.isna().astype(int)
X["f_dosehit"] = X.dose_hit.fillna(0).clip(0, 5)
X["f_tum"] = X.tum_z.fillna(0).clip(-10, 10); X["m_tum"] = X.tum_z.isna().astype(int)
X["f_motif"] = X.motif.fillna(50) / 100
X["f_ks"] = X.ks_rel_pct.fillna(0) / 100           # release (all PSP, in-sample)
X["f_ksdev"] = X.ks_dev_pct.fillna(0) / 100        # DEV-trained, no phospho features
X["test_protein"] = X.test_protein.fillna(0).astype(int)
X["f_ppi"] = X.ppi
NEW = ["f_rep", "f_stouf", "f_ndrop", "f_spec", "f_up", "f_dose", "m_dose", "f_dosehit", "f_tum", "m_tum"]
NEWM = NEW + ["f_motif"]
print(X.shape, X.y_psp.sum(), X.y_ext.sum())


def auc_by_kinase(score, y, kin, min_pos=3):
    out = {}
    for k, idx in pd.Series(np.arange(len(y))).groupby(kin.values).groups.items():
        yy = y[idx]
        if yy.sum() >= min_pos and (1 - yy).sum() > 10:
            out[k] = roc_auc_score(yy, score[idx])
    return out


def metrics(score, y, kin):
    a = auc_by_kinase(score, y, kin)
    top = np.argsort(-score)
    base = y.mean()
    n1 = max(int(0.01 * len(y)), 1)
    return dict(macro_auroc=np.mean(list(a.values())) if a else np.nan, n_kin=len(a), pooled_auroc=roc_auc_score(y, score),
                auprc=average_precision_score(y, score), enrich_top1pct=y[top[:n1]].mean() / base,
                prec_top500=y[top[:500]].mean(), n_pos=int(y.sum()))


def fit_cv(feats, y, groups, C=1.0):
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, groups):
        Z = X[feats].values; mu = Z[tr].mean(0); sd = Z[tr].std(0) + 1e-9
        m = LogisticRegression(C=C, class_weight="balanced", max_iter=2000).fit((Z[tr] - mu) / sd, y[tr])
        oof[te] = m.decision_function((Z[te] - mu) / sd)
    Z = X[feats].values; mu = Z.mean(0); sd = Z.std(0) + 1e-9
    m = LogisticRegression(C=C, class_weight="balanced", max_iter=2000).fit((Z - mu) / sd, y)
    return oof, m.decision_function((Z - mu) / sd), dict(zip(feats, m.coef_[0]))


yP = X.y_psp.values; groups = X.uniprot.values
oof_new, full_new, coef_new = fit_cv(NEW, yP, groups)
oof_newm, full_newm, coef_newm = fit_cv(NEWM, yP, groups)
X["s_new_oof"], X["s_new"] = oof_new, full_new
X["s_newm_oof"], X["s_newm"] = oof_newm, full_newm
INT = NEWM + ["f_ksdev", "f_ppi"]
def fit_dev(feats):
    tr = (X.test_protein == 0).values
    Z = X[feats].values; mu = Z[tr].mean(0); sd = Z[tr].std(0) + 1e-9
    m = LogisticRegression(class_weight="balanced", max_iter=2000).fit((Z[tr] - mu) / sd, yP[tr])
    return m.decision_function((Z - mu) / sd), dict(zip(feats, m.coef_[0]))
X["d_new"], c1 = fit_dev(NEW); X["d_newm"], c2 = fit_dev(NEWM); X["d_int"], c3 = fit_dev(INT)
X["d_ksppi"], c4 = fit_dev(["f_ksdev", "f_ppi"]); X["d_ksmotif"], c5 = fit_dev(["f_ksdev", "f_ppi", "f_motif"])
pd.DataFrame([coef_new, coef_newm, c1, c2, c3, c4, c5], index=["NEW(cv,all)", "NEW+motif(cv,all)", "NEW(dev)", "NEW+motif(dev)", "INT(dev)", "ksdev+ppi(dev)", "ksdev+ppi+motif(dev)"]).T.to_csv(f"{RES}/val_lr_coefficients.tsv", sep="\t")

# ---- PSP evaluation (label-free features; LR OOF)
rows = []
single = {"perturbation: replicated drop (2nd study)": X.f_rep, "perturbation: Stouffer drop": X.f_stouf,
          "perturbation: drop beyond generic": X.f_spec, "dose coherence r": X.f_dose, "tumour co-variation z": X.f_tum,
          "PPI": X.f_ppi, "Cantley motif": X.f_motif, "kinaseSub pct (in-sample on PSP!)": X.f_ks,
          "NEW (LR, OOF)": X.s_new_oof, "NEW+motif (LR, OOF)": X.s_newm_oof}
for nm, s in single.items():
    rows.append(dict(label="PSP", score=nm, **metrics(s.values, yP, X.kinase)))
# subsets with evidence available
for nm, col, msk in [("dose coherence r [dose-evaluable]", "f_dose", X.m_dose == 0), ("tumour z [tumour-evaluable]", "f_tum", X.m_tum == 0)]:
    rows.append(dict(label="PSP", score=nm, **metrics(X.loc[msk, col].values, yP[msk.values], X.kinase[msk])))

T = X[X.test_protein == 1].reset_index(drop=True); yT = T.y_psp.values
for nm, col in [("Cantley motif", "f_motif"), ("kinaseSub DEV-trained, no-phospho", "f_ksdev"), ("PPI", "f_ppi"),
                ("kinaseSub DEV + PPI (LR)", "d_ksppi"), ("kinaseSub DEV + PPI + motif (LR)", "d_ksmotif"),
                ("NEW (LR dev)", "d_new"), ("NEW+motif (LR dev)", "d_newm"), ("INTEGRATED NEW+motif+kinaseSubDEV+PPI (LR dev)", "d_int")]:
    rows.append(dict(label="PSP TEST proteins (held out)", score=nm, **metrics(T[col].values, yT, T.kinase)))
mT = T.ks_dev_pct.notna()
# ---- External evaluation: curated non-PSP edges; PSP positives removed from universe
E = X[X.y_psp == 0].reset_index(drop=True)
yE = E.y_ext.values
# stacked integrated score, grouped OOF on external labels
STK = ["s_newm", "f_ks", "f_ppi"]
oof_int = np.zeros(len(E))
for tr, te in GroupKFold(5).split(E, yE, E.uniprot.values):
    Z = E[STK].values; mu = Z[tr].mean(0); sd = Z[tr].std(0) + 1e-9
    m = LogisticRegression(class_weight="balanced", max_iter=2000).fit((Z[tr] - mu) / sd, yE[tr])
    oof_int[te] = m.decision_function((Z[te] - mu) / sd)
Z = E[STK].values; mu = Z.mean(0); sd = Z.std(0) + 1e-9
mI = LogisticRegression(class_weight="balanced", max_iter=2000).fit((Z - mu) / sd, yE)
print("stack coefs", dict(zip(STK, mI.coef_[0])))
E["s_int_oof"] = oof_int
# simple unfitted rank-sum integration (robustness): mean of within-kinase percentile ranks
for c in ["s_newm", "f_ks", "f_motif", "s_new"]:
    E["r_" + c] = E.groupby("kinase")[c].rank(pct=True)
E["s_rank_int"] = (E.r_s_newm + E.r_f_ks) / 2
ext_scores = {"Cantley motif": E.f_motif, "kinaseSub pct": E.f_ks, "PPI": E.f_ppi,
              "perturbation: replicated drop": E.f_rep, "tumour co-variation z": E.f_tum, "dose coherence r": E.f_dose,
              "NEW (PSP-trained LR)": E.s_new, "NEW+motif (PSP-trained LR)": E.s_newm,
              "integrated: rank-mean(NEW+motif, kinaseSub)": E.s_rank_int, "integrated: stacked (OOF)": E.s_int_oof,
              "kinaseSub DEV-trained, no-phospho": E.f_ksdev, "kinaseSub DEV + PPI + motif (LR)": E.d_ksmotif,
              "INTEGRATED NEW+motif+kinaseSubDEV+PPI (LR dev)": E.d_int}
for nm, s in ext_scores.items():
    rows.append(dict(label="external curated non-PSP", score=nm, **metrics(s.values, yE, E.kinase)))
# also text-mining tier
for nm, s in [("kinaseSub pct", E.f_ks), ("Cantley motif", E.f_motif), ("NEW+motif (PSP-trained LR)", E.s_newm),
              ("integrated: rank-mean(NEW+motif, kinaseSub)", E.s_rank_int), ("integrated: stacked (OOF)", E.s_int_oof)]:
    msk = (E.y_ext == 0).values
    rows.append(dict(label="external text-mining non-PSP", score=nm, **metrics(s.values[msk], E.y_ext_txt.values[msk], E.kinase[msk])))
V = pd.DataFrame(rows); V.to_csv(f"{RES}/val_metrics.tsv", sep="\t", index=False)
print(V.round(3).to_string())

# ---- bootstrap CI (by kinase) for key external comparison
def boot_delta(a, b, y, kin, B=300):
    ak = auc_by_kinase(a, y, kin); bk = auc_by_kinase(b, y, kin)
    ks = sorted(set(ak) & set(bk)); d = np.array([ak[k] - bk[k] for k in ks])
    bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(B)]
    return d.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(ks)
brow = []
for a, b in [("d_int", "d_ksmotif"), ("d_int", "f_ksdev"), ("d_int", "f_motif"), ("d_newm", "f_motif")]:
    m, lo, hi, n = boot_delta(T[a].values, T[b].values, yT, T.kinase)
    brow.append(dict(comparison=f"{a} - {b}", label="PSP TEST", d_macro_auroc=m, lo=lo, hi=hi, n_kin=n))
    m, lo, hi, n = boot_delta(E[a].values, E[b].values, yE, E.kinase)
    brow.append(dict(comparison=f"{a} - {b}", label="external", d_macro_auroc=m, lo=lo, hi=hi, n_kin=n))
for a, b in [("integrated: rank-mean(NEW+motif, kinaseSub)", "kinaseSub pct"), ("integrated: stacked (OOF)", "kinaseSub pct"),
             ("NEW+motif (PSP-trained LR)", "Cantley motif"), ("integrated: rank-mean(NEW+motif, kinaseSub)", "Cantley motif")]:
    m, lo, hi, n = boot_delta(ext_scores[a].values, ext_scores[b].values, yE, E.kinase)
    brow.append(dict(comparison=f"{a} - {b}", label="external", d_macro_auroc=m, lo=lo, hi=hi, n_kin=n))
m, lo, hi, n = boot_delta(X.s_newm_oof.values, X.f_motif.values, yP, X.kinase)
brow.append(dict(comparison="NEW+motif (OOF) - Cantley motif", label="PSP", d_macro_auroc=m, lo=lo, hi=hi, n_kin=n))
m, lo, hi, n = boot_delta(X.s_new_oof.values, -X.f_ppi.values * 0 + rng.random(len(X)), yP, X.kinase)
brow.append(dict(comparison="NEW (OOF) - random", label="PSP", d_macro_auroc=m, lo=lo, hi=hi, n_kin=n))
B = pd.DataFrame(brow); B.to_csv(f"{RES}/val_bootstrap.tsv", sep="\t", index=False); print(B.round(3).to_string())

# ---- independence of new evidence from kinaseSub
ind = []
for c in ["f_rep", "f_stouf", "f_tum", "f_dose", "s_new"]:
    m = X.ks_rel_pct.notna()
    ind.append(dict(feature=c, spearman_vs_ks=X.loc[m, c].corr(X.loc[m, "ks_rel_pct"], method="spearman"),
                    spearman_vs_motif=X[c].corr(X.f_motif, method="spearman")))
# NEW evidence AUROC within kinaseSub strata (does it add information beyond kinaseSub?)
X["ks_bin"] = pd.cut(X.ks_rel_pct.fillna(0), [-1, 50, 90, 99, 100.1], labels=["absent/<50", "50-90", "90-99", ">99"])
for b, d in X.groupby("ks_bin", observed=True):
    if d.y_psp.sum() >= 20:
        ind.append(dict(feature=f"NEW OOF AUROC within kinaseSub bin {b}", spearman_vs_ks=roc_auc_score(d.y_psp, d.s_new_oof),
                        spearman_vs_motif=d.y_psp.sum()))
Eb = E.copy(); Eb["ks_bin"] = pd.cut(Eb.ks_rel_pct.fillna(0), [-1, 50, 90, 99, 100.1], labels=["absent/<50", "50-90", "90-99", ">99"])
for b, d in Eb.groupby("ks_bin", observed=True):
    if d.y_ext.sum() >= 10:
        ind.append(dict(feature=f"[external] NEW AUROC within kinaseSub bin {b}", spearman_vs_ks=roc_auc_score(d.y_ext, d.s_new),
                        spearman_vs_motif=d.y_ext.sum()))
I = pd.DataFrame(ind); I.to_csv(f"{RES}/val_independence.tsv", sep="\t", index=False); print(I.round(3).to_string())

# ---- site-centric: kinase-level vs family-level accuracy
def site_centric(df, score, ycol):
    pos_sites = df.loc[df[ycol] == 1, "site"].unique()
    d = df[df.site.isin(pos_sites)][["site", "kinase", "family", ycol, score]]
    d = d[d.groupby("site").kinase.transform("nunique") >= 3]
    top = d.sort_values(score, ascending=False).groupby("site").head(1)
    fam_true = d[d[ycol] == 1].groupby("site").family.apply(set)
    exact = top[ycol].mean()
    fam = np.mean([f in fam_true.get(s, set()) for s, f in zip(top.site, top.family)])
    # chance: random kinase among candidates
    ch_exact = d.groupby("site")[ycol].mean().mean()
    ch_fam = np.mean([np.mean([f in fam_true[s] for f in g.family]) for s, g in d.groupby("site")])
    return dict(n_sites=len(top), top1_exact=exact, top1_family=fam, chance_exact=ch_exact, chance_family=ch_fam)
sc = []
for nm, col in [("NEW (OOF)", "s_new_oof"), ("NEW+motif (OOF)", "s_newm_oof"), ("Cantley motif", "f_motif"), ("kinaseSub (in-sample)", "f_ks")]:
    sc.append(dict(label="PSP", score=nm, **site_centric(X, col, "y_psp")))
for nm, col in [("d_int", "d_int"), ("d_ksmotif", "d_ksmotif"), ("f_ksdev", "f_ksdev")]:
    sc.append(dict(label="PSP TEST", score=nm, **site_centric(T, col, "y_psp")))
for nm, col in [("NEW", "s_new"), ("NEW+motif", "s_newm"), ("Cantley motif", "f_motif"), ("kinaseSub", "f_ks"), ("integrated stacked OOF", "s_int_oof"),
                ("d_int", "d_int"), ("d_ksmotif", "d_ksmotif"), ("f_ksdev", "f_ksdev"),
                ("integrated rank-mean", "s_rank_int")]:
    sc.append(dict(label="external", score=nm, **site_centric(E, col, "y_ext")))
SC = pd.DataFrame(sc); SC.to_csv(f"{RES}/val_site_centric_family.tsv", sep="\t", index=False); print(SC.round(3).to_string())

# ---- kinase-label permutation and family-restricted swaps (PSP, NEW OOF score)
fam = X.drop_duplicates("kinase").set_index("kinase").family
kins = np.array(sorted(X.kinase.unique()))
psp_pairs = set(zip(*[X.kinase[X.y_psp == 1], X.site[X.y_psp == 1]]))
PSPall = pd.read_parquet(f"{OUT}/psp_edges.parquet"); psp_all = set(zip(PSPall.kinase, PSPall.site))
def perm_auc(mapping, score="s_new_oof"):
    yy = np.array([(mapping[k], s) in psp_all for k, s in zip(X.kinase, X.site)], dtype=int)
    a = auc_by_kinase(X[score].values, yy, X.kinase)
    return np.mean(list(a.values())), len(a)
pr = [dict(control="observed", score=s, macro_auroc=perm_auc({k: k for k in kins}, s)[0]) for s in ["s_new_oof", "s_newm_oof"]]
for i in range(20):
    p = dict(zip(kins, rng.permutation(kins)))
    for s in ["s_new_oof", "s_newm_oof"]:
        pr.append(dict(control="random kinase permutation", score=s, macro_auroc=perm_auc(p, s)[0]))
    # swap within family (derangement where possible), across family
    pw = {}
    for f, ks in pd.Series(kins).groupby(fam.reindex(kins).values):
        ks = list(ks)
        sh = list(np.roll(rng.permutation(ks), 1)) if len(ks) > 1 else ks
        pw.update(dict(zip(ks, sh)))
    moved = [k for k in kins if pw[k] != k]
    for s in ["s_new_oof", "s_newm_oof"]:
        yy = np.array([(pw[k], si) in psp_all for k, si in zip(X.kinase, X.site)], dtype=int)
        msk = X.kinase.isin(moved).values
        a = auc_by_kinase(X[s].values[msk], yy[msk], X.kinase[msk])
        a0 = auc_by_kinase(X[s].values[msk], X.y_psp.values[msk], X.kinase[msk])
        pr.append(dict(control="within-family swap (moved kinases)", score=s, macro_auroc=np.mean(list(a.values())),
                       observed_same_kinases=np.mean(list(a0.values()))))
PR = pd.DataFrame(pr)
PRs = PR.groupby(["control", "score"]).agg(mean=("macro_auroc", "mean"), sd=("macro_auroc", "std"), observed_same=("observed_same_kinases", "mean")).reset_index()
PRs.to_csv(f"{RES}/val_permutation_controls.tsv", sep="\t", index=False); print(PRs.round(3).to_string())

# per-kinase AUROC table
pk = []
for nm, col in [("NEW_oof", "s_new_oof"), ("NEWmotif_oof", "s_newm_oof"), ("motif", "f_motif"), ("rep_drop", "f_rep"), ("tumour", "f_tum")]:
    for k, v in auc_by_kinase(X[col].values, yP, X.kinase).items():
        pk.append(dict(kinase=k, score=nm, auroc=v))
PK = pd.DataFrame(pk).pivot(index="kinase", columns="score", values="auroc")
PK["n_psp"] = X.groupby("kinase").y_psp.sum().reindex(PK.index)
PK.to_csv(f"{RES}/val_per_kinase_psp.tsv", sep="\t"); print(PK.round(2).to_string())

X.to_parquet(f"{OUT}/scored.parquet"); E[["kinase", "site", "s_int_oof", "s_rank_int"]].to_parquet(f"{OUT}/scored_ext.parquet")
