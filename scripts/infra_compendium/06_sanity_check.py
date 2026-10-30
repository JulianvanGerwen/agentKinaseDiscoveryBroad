"""Quick sanity check: in kinase-inhibitor conditions, are the known targets' activities decreased?

Conditions: perturbation_class == kinase_inhibitor, subclass in {inhibitor, inhibitor+stimulus, inhibitor_reversed_contrast
(z sign flipped)}, target_confidence in {high, medium}, qc_flag == ok, dose_nM >= 100 or unknown.
Targets: inhibitor_primary_targets_in_universe (curated primary targets).
Per condition & method: AUROC (targets ranked as most decreased vs all other scored kinases), fraction of target z < 0.
Output: results/infra_compendium/sanity_inhibitor_targets.tsv (per condition x method) and _summary.tsv
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from common import RES

A = pd.read_csv(RES / "conditions_annotated.tsv", sep="\t", dtype={"study_pmid": str})
sel = A[A.perturbation_class.eq("kinase_inhibitor")
        & A.perturbation_subclass.isin(["inhibitor", "inhibitor+stimulus", "inhibitor_reversed_contrast"])
        & A.target_confidence.isin(["high", "medium"]) & A.qc_flag.eq("ok")
        & (A.dose_nM.isna() | A.dose_nM.ge(100))
        & A.inhibitor_primary_targets_in_universe.notna()]
print(f"{len(sel)} inhibitor conditions from {sel.study_pmid.nunique()} studies")
L = pd.read_parquet(RES / "activity_long.parquet", columns=["condition", "kinase", "method", "z", "n_substrates", "z_core"])
L = L[L.condition.isin(sel.conditionID)]
sign = dict(zip(sel.conditionID, np.where(sel.perturbation_subclass.eq("inhibitor_reversed_contrast"), -1.0, 1.0)))
tg = dict(zip(sel.conditionID, sel.inhibitor_primary_targets_in_universe.str.split(";")))
rows = []
variants = {"corev3": lambda d: d[d.method == "corev3"],
            "corev3_nsub>=1": lambda d: d[(d.method == "corev3") & (d.n_substrates >= 1)],
            "corev3_zcore_only": lambda d: d[d.method == "corev3"].assign(z=lambda q: q.z_core),
            "psp_z": lambda d: d[d.method == "psp_z"]}
for (c), d in L.groupby("condition"):
    for name, f in variants.items():
        x = f(d)
        if x.empty:
            continue
        z = x.z.to_numpy() * sign[c]
        is_t = x.kinase.isin(tg[c]).to_numpy()
        nt, nn = is_t.sum(), (~is_t).sum()
        if nt == 0 or nn == 0:
            continue
        r = rankdata(-z)                         # most decreased gets the highest rank
        auc = (r[is_t].sum() - nt * (nt + 1) / 2) / (nt * nn)
        rows.append(dict(condition=c, method=name, n_targets_scored=int(nt), n_kinases=int(len(z)), auroc=auc,
                         frac_targets_negative=float((z[is_t] < 0).mean()), mean_target_z=float(z[is_t].mean()),
                         n_targets_z_lt_m1_96=int((z[is_t] < -1.96).sum())))
R = pd.DataFrame(rows).merge(sel[["conditionID", "study_pmid", "kinase_inhibitors", "inhibitor_primary_targets_in_universe",
                                  "target_confidence"]], left_on="condition", right_on="conditionID").drop(columns="conditionID")
R.to_csv(RES / "sanity_inhibitor_targets.tsv", sep="\t", index=False, float_format="%.4g")


def summ(g):
    pooled_neg = (g.frac_targets_negative * g.n_targets_scored).sum() / g.n_targets_scored.sum()
    return pd.Series(dict(n_conditions=len(g), n_studies=g.study_pmid.nunique(), median_auroc=g.auroc.median(),
                          mean_auroc=g.auroc.mean(), study_balanced_median_auroc=g.groupby("study_pmid").auroc.median().median(), frac_cond_auroc_gt_0_5=(g.auroc > 0.5).mean(),
                          pooled_frac_target_z_negative=pooled_neg, median_target_z=g.mean_target_z.median()))


S = R.groupby("method").apply(summ)
# same condition set across methods
both = set(R[R.method == "psp_z"].condition) & set(R[R.method == "corev3"].condition)
S2 = R[R.condition.isin(both)].groupby("method").apply(summ)
S2.index = [f"{i} (shared conds)" for i in S2.index]
S3 = R[R.target_confidence.eq("high")].groupby("method").apply(summ)
S3.index = [f"{i} (high-conf targets only)" for i in S3.index]
out = pd.concat([S, S2, S3])
out.to_csv(RES / "sanity_inhibitor_targets_summary.tsv", sep="\t", float_format="%.3f")
print(out.round(3).to_string())
