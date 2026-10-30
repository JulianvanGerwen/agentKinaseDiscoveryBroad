"""I03 step 4b: full kinaseSub scores for panel kinases (all sites, not only top200/top5).
 - ks_rel_pct : release model (trained on DEV+TEST = all PSP; includes F6 phospho features) -> in-sample on PSP.
 - ks_dev_pct : 'nophospho' DEV-trained R2 (no compendium/CPTAC features, TEST proteins never seen) -> fair on TEST PSP
                and genuinely independent of the perturbation/tumour evidence.
 - test_protein: substrate protein in kinaseSub TEST split.
Output cache/b_I03/ks_full.parquet
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow.dataset as ds
KSM = "<HACKATHON_ROOT>/kinaseSubMain"
OUT = "/workspace/cache/b_I03"
KC = pd.read_csv(f"{OUT}/kinase_conditions.tsv", sep="\t")
panel = sorted(set(KC.kinase) - {"PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"})
kin = pd.read_parquet(f"{KSM}/cache/index/kinases.parquet")
sites = pd.read_parquet(f"{KSM}/cache/index/sites.parquet", columns=["site_id", "acc"])
kidx = {g: i for i, g in enumerate(kin.gene) if g in panel}
print(len(kidx), "panel kinases in kinaseSub")
test = set(pd.read_csv(f"{KSM}/results/splits/test_proteins.tsv", sep="\t").acc)
comp = set(pd.read_parquet(f"{OUT}/site_annot.parquet", columns=["site"]).site)

def load(path, col):
    d = ds.dataset(path)
    t = d.to_table(columns=["kinase_idx", "site_idx", col], filter=pc.field("kinase_idx").isin(list(kidx.values())))
    df = t.to_pandas()
    df["pct"] = df.groupby("kinase_idx")[col].rank(pct=True) * 100
    return df
a = load(f"{KSM}/cache/release/r2r3_release__20261027-012616__s0/release_scores.parquet", "pct").rename(columns={"pct": "ks_rel_pct"})
a = a[["kinase_idx", "site_idx", "ks_rel_pct"]] if "ks_rel_pct" in a else a
b = load(f"{KSM}/cache/runs/r2_final_nophospho_uobs__20261026-213950__s0/uobs_predictions.parquet", "score")
b = b.rename(columns={"pct": "ks_dev_pct", "score": "ks_dev_score"})
a = a.loc[:, ~a.columns.duplicated()]
M = a.merge(b, on=["kinase_idx", "site_idx"], how="outer")
inv = {i: g for g, i in kidx.items()}
M["kinase"] = M.kinase_idx.map(inv)
M["site"] = sites.site_id.values[M.site_idx.values]
M = M[M.site.isin(comp)]
M["test_protein"] = sites.acc.values[M.site_idx.values]
M["test_protein"] = M.test_protein.isin(test).astype(int)
M = M[["kinase", "site", "ks_rel_pct", "ks_dev_pct", "ks_dev_score", "test_protein"]]
M.to_parquet(f"{OUT}/ks_full.parquet"); print(M.shape, M.isna().mean().to_dict())
