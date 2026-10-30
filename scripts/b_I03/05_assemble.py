"""I03 step 5: assemble the kinase x site evidence table (panel kinases, residue-class matched, perturbation-evaluable).
Output cache/b_I03/pairs.parquet
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
OUT = "/workspace/cache/b_I03"
P = pd.read_parquet(f"{OUT}/pert_long.parquet"); P["kinase"] = P.kinase.astype(str)
LIPID = {"PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"}
P = P[~P.kinase.isin(LIPID)]
KA = pd.read_csv("/workspace/results/infra_compendium/kinase_annotation.tsv", sep="\t").set_index("kinase")
SA = pd.read_parquet(f"{OUT}/site_annot.parquet").set_index("site")
P["residue"] = SA.residue.reindex(P.site).values
rc = KA.residue_class.reindex(P.kinase).fillna("ST").values
okres = np.where(pd.Series(rc).str.contains("Y").values & pd.Series(rc).str.contains("S").values, True,
                 np.where(pd.Series(rc).str.contains("Y").values, P.residue.values == "Y", P.residue.values != "Y"))
P = P[okres]
print(KA.residue_class.value_counts().to_dict())
D = pd.read_parquet(f"{OUT}/dose_long.parquet")
T = pd.read_parquet(f"{OUT}/tumour_long.parquet")
KS = pd.read_parquet(f"{OUT}/ks_pairs.parquet")
KF = pd.read_parquet(f"{OUT}/ks_full.parquet")
X = P.merge(D, on=["kinase", "site"], how="left").merge(T, on=["kinase", "site"], how="left").merge(KS, on=["kinase", "site"], how="left").merge(KF, on=["kinase", "site"], how="left")
# motif
mot = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
mcols = set(mot.columns)
X["motif"] = np.nan
for k, idx in X.groupby("kinase").groups.items():
    if k in mcols:
        X.loc[idx, "motif"] = mot[k].reindex(X.loc[idx, "site"]).values
# PPI (kinase gene - substrate gene)
pp = pd.read_parquet(f"{OUT}/ppi_pairs.parquet")
pps = set(zip(pp.a, pp.b))
X["sub_gene"] = SA.gene.reindex(X.site).astype(str).str.upper().values
X["ppi"] = [int((k, g) in pps) for k, g in zip(X.kinase, X.sub_gene)]
# labels
PSP = pd.read_parquet(f"{OUT}/psp_edges.parquet"); EX = pd.read_parquet(f"{OUT}/ext_edges.parquet"); KN = pd.read_parquet(f"{OUT}/known_any.parquet")
key = X.kinase + "|" + X.site
X["y_psp"] = key.isin(set(PSP.kinase + "|" + PSP.site)).astype(int)
exc = EX[EX.tier == "curated"]
X["y_ext"] = key.isin(set(exc.kinase + "|" + exc.site)).astype(int)
X["y_ext_txt"] = key.isin(set(EX[EX.tier == "textmining"].kinase + "|" + EX[EX.tier == "textmining"].site)).astype(int)
X["known_any"] = key.isin(set(KN.kinase + "|" + KN.site)).astype(int)
# any PSP kinase at this site (site has a known kinase)
X["site_has_psp"] = X.site.isin(set(PSP.site)).astype(int)
X["family"] = KA.manning_family.reindex(X.kinase).values
X["group"] = KA.manning_group.reindex(X.kinase).values
for c in ["functional_score", "psp_regulatory", "psp_reg_function", "gene", "n_studies"]:
    X[c] = SA[c].reindex(X.site).values
X.to_parquet(f"{OUT}/pairs.parquet")
print(X.shape, X.kinase.nunique())
print(X.groupby("kinase").agg(n=("site", "size"), rep=("n_st", lambda x: (x >= 2).sum()), psp=("y_psp", "sum"), ext=("y_ext", "sum"),
                              dose=("dose_r", lambda x: x.notna().sum()), tum=("tum_z", lambda x: x.notna().sum())).to_string())
