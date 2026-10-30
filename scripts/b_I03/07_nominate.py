"""I03 step 7: nominate new kinase->site edges.
Score (log-odds style, weights from the conditional model 06c, PSP all proteins):
  prior  = logit of PSP rate in the pair's kinaseSub-DEV tier (from val_stratified table, no-drop rows) -- sequence/PPI prior
  + log OR(drop in >=2 studies) [+ log OR(drop in 1 study)] + log OR(tumour z>3) + log OR(dose r>0.3)
Hard filters: not known in PSP/OmniPath(any source)/iPTMnet; decrease in >=2 independent studies (PMIDs) of K-loss;
no study with an increase; drop beyond the site's generic inhibitor response (f_spec>0); kinaseSub-DEV pct>=95 or motif>=90.
Controls per edge:
  rank_inh   : rank of K among all panel kinases with >=2 studies at this site by Stouffer drop (1 = K's inhibitors give the
               strongest drop); kinases sharing the identical inhibitor set are tied (reported as ties)
  best_other_family : best-scoring family member and its kinaseSub-DEV pct (family-swap control)
  p_cov      : coverage-matched random-site p: fraction of the kinase's sites measured in the same number of studies whose
               Stouffer drop is at least as strong
Outputs results/b_I03/tables/edges_ranked_all.tsv(.gz), edges_novel_candidates.tsv
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
OUT = "/workspace/cache/b_I03"; RES = "/workspace/results/b_I03/tables"
X = pd.read_parquet(f"{OUT}/scored.parquet")
G = pd.read_parquet(f"{OUT}/gbm_scores.parquet")
X = X.merge(G[["kinase", "site", "g_PRIOR", "g_FULL"]], on=["kinase", "site"], how="left")
OR = pd.read_csv(f"{RES}/val_conditional_OR.tsv", sep="\t")
OR = OR[OR.label == "PSP, all proteins"].set_index("term").OR
ST = pd.read_csv(f"{RES}/val_stratified_psp_rate.tsv", sep="\t")
pri = ST[ST.pert == "no drop"].groupby("ks_tier").apply(lambda d: (d.psp_rate * d.n).sum() / d.n.sum(), include_groups=False)
X["ks_tier"] = pd.cut(X.ks_dev_pct.fillna(0), [-1, 95, 99, 99.7, 100.1], labels=["<95", "95-99", "99-99.7", ">99.7"]).astype(str)
p0 = X.ks_tier.map(pri).astype(float).clip(1e-4, 0.5)
X["prior_logit"] = np.log(p0 / (1 - p0))
X["drop2"] = X.n_st_drop >= 2; X["drop1"] = X.n_st_drop == 1
X["tum_pos"] = X.tum_z > 3; X["dose_pos"] = X.dose_r > 0.3
X["score"] = (X.prior_logit + np.log(OR["drop2"]) * X.drop2 + np.log(OR["drop1"]) * X.drop1
              + np.log(OR["tum_pos"]) * X.tum_pos + np.log(OR["dose_pos"]) * X.dose_pos)
X["post_prob"] = 1 / (1 + np.exp(-X.score))

# rank of K's inhibitor drop among panel kinases at the site
X["rank_inh"] = X.groupby("site").stouffer.rank(method="min")
X["n_kin_site"] = X.groupby("site").kinase.transform("size")
# coverage-matched random-site empirical p (same kinase, same n_st)
X["p_cov"] = X.groupby(["kinase", "n_st"]).stouffer.rank(method="max", pct=True)
# best other family member by kinaseSub-DEV (family swap)
fam_best = []
for (s, f), d in X.groupby(["site", "family"]):
    if len(d) > 1:
        for r in d.itertuples():
            o = d[d.kinase != r.kinase].sort_values("ks_dev_pct", ascending=False).iloc[0]
            fam_best.append((r.kinase, s, o.kinase, o.ks_dev_pct, o.stouffer))
FB = pd.DataFrame(fam_best, columns=["kinase", "site", "fam_best_other", "fam_other_ksdev", "fam_other_stouffer"])
X = X.merge(FB, on=["kinase", "site"], how="left")
cols = ["kinase", "family", "site", "gene", "sequence_window", "score", "post_prob", "g_FULL", "ks_dev_pct", "ks_rel_pct", "motif", "ppi",
        "n_st", "n_st_drop", "n_st_up", "n_cond", "stouffer", "mean_z", "min2_z", "gen_z", "p_cov", "rank_inh", "n_kin_site",
        "dose_r", "dose_hit", "dose_nstrong", "dose_gap", "tum_z", "tum_ncoh", "tum_npos", "tum_rmed",
        "fam_best_other", "fam_other_ksdev", "fam_other_stouffer", "y_psp", "y_ext", "y_ext_txt", "known_any", "site_has_psp",
        "psp_regulatory", "psp_reg_function", "functional_score", "test_protein"]
SA = pd.read_parquet(f"{OUT}/site_annot.parquet").set_index("site")
X["sequence_window"] = SA.sequence_window.reindex(X.site).values
X = X.sort_values("score", ascending=False)
X[cols].to_csv(f"{RES}/edges_ranked_all.tsv.gz", sep="\t", index=False, float_format="%.4g")
# family-level novelty: no kinase of the same Manning family known for this site in PSP/OmniPath/iPTMnet
KN = pd.read_parquet(f"{OUT}/known_any.parquet")
KA = pd.read_csv("/workspace/results/infra_compendium/kinase_annotation.tsv", sep="\t").set_index("kinase")
KN["family"] = KA.manning_family.reindex(KN.kinase).values
fam_known = set(zip(KN.family.dropna(), KN.site[KN.family.notna()]))
site_known = KN.groupby("site").kinase.apply(lambda x: ";".join(sorted(set(x))))
X["family_known"] = [(f, s) in fam_known for f, s in zip(X.family, X.site)]
X["site_known_kinases"] = site_known.reindex(X.site).values
X["orphan_site"] = X.site_known_kinases.isna()
# explaining-away: among panel kinases whose loss reduces this site (>=1 study), is K the best sequence/PPI-supported one?
X["drop_any"] = X.n_st_drop >= 1
X["explain_rank"] = np.nan
m = X.drop_any
X.loc[m, "explain_rank"] = X[m].groupby("site").ks_dev_pct.rank(ascending=False, method="min")
X["n_drop_kinases"] = X.groupby("site").drop_any.transform("sum")
X["self_site"] = X.gene.astype(str).str.upper() == X.kinase
EXCL_K = {"MAP2K1", "MAP2K2", "BRAF", "RAF1", "ARAF"}   # MEK/RAF: single-substrate kinases -> any 'new' substrate is downstream
cols += ["family_known", "orphan_site", "site_known_kinases", "explain_rank", "n_drop_kinases"]
X[cols].to_csv(f"{RES}/edges_ranked_all.tsv.gz", sep="\t", index=False, float_format="%.4g")
C = X[(~X.family_known) & (X.known_any == 0) & (X.n_st_drop >= 2) & (X.n_st_up == 0) & (X.f_spec > 0) & ((X.ks_dev_pct >= 95) | (X.motif >= 90))
      & (X.explain_rank <= 1) & (~X.self_site) & (~X.kinase.isin(EXCL_K))]
C = C.sort_values(["score", "stouffer"], ascending=[False, True])
C[cols].to_csv(f"{RES}/edges_novel_candidates.tsv", sep="\t", index=False, float_format="%.4g")
print(len(X), "pairs;", len(C), "novel candidates")
print(C.groupby("kinase").size().sort_values(ascending=False).head(30).to_dict())
show = ["kinase", "site", "gene", "orphan_site", "n_drop_kinases", "score", "ks_dev_pct", "motif", "n_st", "n_st_drop", "stouffer", "p_cov", "rank_inh", "dose_r", "tum_z",
        "fam_best_other", "fam_other_ksdev", "psp_regulatory", "functional_score"]
pd.set_option("display.width", 250)
print(C[show].head(120).round(2).to_string())
