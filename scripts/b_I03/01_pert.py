"""I03 step 1: perturbation (primary-direction) evidence per kinase x site.
For each kinase K: 'K-loss' conditions = kinase-inhibitor conditions with K among curated primary targets
(qc ok, single agent, not inhibitor-vs-inhibitor, not pan/antibody/indirect agents, dose >=100 nM or unknown,
<=6 primary targets) + genetic K(-) conditions. Site log2FC -> robust z per condition (median/MAD); sign flipped
for reversed contrasts. Per study (PMID): mean z over that study's K-loss conditions. Only DECREASES are counted
as evidence (reverse-direction responses are not portable across labs).
Outputs cache/b_I03/pert_long.parquet (kinase, site, n_st, n_st_drop, stouffer_neg, mean_z, min2_z, gen_z, ...)
and cache/b_I03/kinase_conditions.tsv.
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/cache/b_I03"

C = load.conditions()
EXCL = {"SHP099", "staurosporine", "pertuzumab", "trastuzumab", "perifosine"}
inh = C[(C.perturbation_class == "kinase_inhibitor") & (C.qc_flag == "ok") & (~C.is_combination.astype(bool))
        & (C.perturbation_subclass != "inhibitor_vs_other_inhibitor") & C.inhibitor_primary_targets.notna()
        & (~C.agents.isin(EXCL))].copy()
inh = inh[(inh.dose_nM.isna()) | (inh.dose_nM >= 100)]
inh["targets"] = inh.inhibitor_primary_targets.str.split(";")
inh = inh[inh.targets.str.len() <= 6]
rows = []
for r in inh.itertuples():
    for t in r.targets:
        rows.append((t, r.conditionID, r.study_pmid, r.agents, -1.0 if r.contrast_reversed else 1.0, "inhibitor"))
gen = C[C.genetic_target_kinases.notna() & (C.qc_flag != "scale_suspect")]
for r in gen.itertuples():
    for t in str(r.genetic_target_kinases).split(";"):
        if t.endswith("(-)"):
            rows.append((t[:-3], r.conditionID, r.study_pmid, "genetic:" + t, -1.0 if r.contrast_reversed else 1.0, "genetic"))
KC = pd.DataFrame(rows, columns=["kinase", "conditionID", "study", "agent", "sign", "kind"]).drop_duplicates(["kinase", "conditionID"])
ns = KC.groupby("kinase").study.nunique()
panel = sorted(ns[ns >= 2].index)
KC = KC[KC.kinase.isin(panel)]
KC.to_csv(f"{OUT}/kinase_conditions.tsv", sep="\t", index=False)
print(len(panel), "kinases with >=2 studies:", panel)

# all inhibitor conditions (for the generic-decrease control)
allinh = sorted(set(inh.conditionID))
need = sorted(set(KC.conditionID) | set(allinh))
X = load.sites_matrix(need)
med = X.median(axis=0); mad = (X - med).abs().median(axis=0) * 1.4826
Z = ((X - med) / mad.clip(lower=0.05)).astype(np.float32)
sign = KC.drop_duplicates("conditionID").set_index("conditionID").sign
for c in sign.index[sign < 0]:
    Z[c] = -Z[c]
Z = Z.clip(-8, 8)
sites = Z.index.values
print("Z", Z.shape)

# generic decrease tendency per site: mean over studies of mean z across all inhibitor conditions
cst = C.set_index("conditionID").study_pmid
def study_mean(cols):
    st = cst.reindex(cols).values
    out = {}
    for s in np.unique(st):
        cc = [c for c, x in zip(cols, st) if x == s]
        out[s] = Z[cc].mean(axis=1, skipna=True)
    return pd.DataFrame(out)

out = []
for k in panel:
    kc = KC[KC.kinase == k]
    M = study_mean(list(kc.conditionID))            # sites x studies
    # generic control: all inhibitor conditions whose targets do NOT include k
    other = [c for c in allinh if c not in set(kc.conditionID)]
    G = study_mean(other)
    gen_z = G.mean(axis=1, skipna=True)
    n_st = M.notna().sum(axis=1)
    keep = n_st >= 1
    Mv = M[keep]
    arr = Mv.values
    neg = np.where(np.isnan(arr), np.nan, np.minimum(arr, 0))
    n = (~np.isnan(arr)).sum(1)
    stouf = np.nansum(arr, 1) / np.sqrt(n)
    srt = np.sort(np.where(np.isnan(arr), np.inf, arr), 1)
    min2 = np.where(n >= 2, srt[:, 1], np.nan)       # 2nd most negative study mean (replication)
    df = pd.DataFrame(dict(kinase=k, site=Mv.index.values, n_st=n,
                           n_st_drop=(arr < -1.0).sum(1), n_st_up=(arr > 1.0).sum(1),
                           stouffer=stouf.astype(np.float32), mean_z=np.nanmean(arr, 1).astype(np.float32),
                           min1_z=srt[:, 0].astype(np.float32), min2_z=min2.astype(np.float32),
                           gen_z=gen_z.reindex(Mv.index).values.astype(np.float32),
                           n_cond=Z.loc[Mv.index, kc.conditionID].notna().sum(1).values))
    out.append(df)
    print(k, len(kc), kc.study.nunique(), len(df), (df.n_st >= 2).sum())
P = pd.concat(out, ignore_index=True)
P["kinase"] = P.kinase.astype("category")
P.to_parquet(f"{OUT}/pert_long.parquet")
print(P.shape)
