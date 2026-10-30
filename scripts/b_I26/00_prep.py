"""I26 step 0: inhibitor-resolvable kinase groups, LoF conditions, annotated PSP edge library, robust z matrix.
Reuses I03 LoF condition definition (cache/b_I03/kinase_conditions.tsv: curated primary target of a single-agent
inhibitor, qc ok, dose>=100 nM/unknown, <=6 primary targets; or genetic K(-)).
Outputs cache/b_I26/: groups.tsv, lof_conditions.tsv, psp_edges_annot.parquet, Z_lof.parquet (sites x conds, robust z,
sign-flipped for reversed contrasts, clipped +-8).
"""
import os, sys, re
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
from common import read_psp
OUT = "/workspace/cache/b_I26"
D = "/workspace/data"

GROUPS = {
    "AKT": ["AKT1", "AKT2", "AKT3"], "ERK": ["MAPK1", "MAPK3"], "MEK": ["MAP2K1", "MAP2K2"],
    "RAF": ["BRAF", "RAF1", "ARAF"], "P38": ["MAPK11", "MAPK14"], "JNK": ["MAPK8", "MAPK9", "MAPK10"],
    "CK2": ["CSNK2A1", "CSNK2A2"], "CAMK2": ["CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G"],
    "PKC": ["PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCZ"], "GSK3": ["GSK3A", "GSK3B"],
    "CDK1/2/5/9": ["CDK1", "CDK2", "CDK5", "CDK9"], "CDK4/6": ["CDK4", "CDK6"], "CDK7/12/13": ["CDK7", "CDK12", "CDK13"],
    "JAK": ["JAK1", "JAK2", "JAK3"], "TBK1/IKKE": ["TBK1", "IKBKE"], "ROCK": ["ROCK1", "ROCK2"],
    "ABL/KIT/PDGFR": ["ABL1", "KIT", "PDGFRA", "PDGFRB"], "ALK/ROS1": ["ALK", "ROS1"], "MET": ["MET"],
    "EGFR/ERBB2": ["EGFR", "ERBB2"], "MTOR": ["MTOR"], "S6K": ["RPS6KB1"], "ATM": ["ATM"], "ATR": ["ATR"],
    "DNAPK": ["PRKDC"], "PLK1": ["PLK1"], "IKKB": ["IKBKB"], "TAK1": ["MAP3K7"], "FLT3": ["FLT3"],
}
k2g = {k: g for g, ks in GROUPS.items() for k in ks}
pd.DataFrame([(g, k) for g, ks in GROUPS.items() for k in ks], columns=["group", "kinase"]).to_csv(f"{OUT}/groups.tsv", sep="\t", index=False)

KC = pd.read_csv("/workspace/cache/b_I03/kinase_conditions.tsv", sep="\t", dtype={"study": str})
KC = KC[KC.kinase.isin(k2g)].assign(group=lambda d: d.kinase.map(k2g))
GC = KC.drop_duplicates(["group", "conditionID"])[["group", "conditionID", "study", "agent", "sign", "kind"]]
GC.to_csv(f"{OUT}/lof_conditions.tsv", sep="\t", index=False)
print(GC.groupby("group").agg(n_cond=("conditionID", "size"), n_st=("study", "nunique")).to_string())

# ---------------- PSP edge library (same as infra psp_z library) + annotation
lib = pd.read_csv("/workspace/results/infra_compendium/psp_ks_edges_human_mapped.tsv", sep="\t")
ks = read_psp(f"{D}/Kinase_Substrate_Dataset")
ps = read_psp(f"{D}/Phosphorylation_site_dataset")
psh = ps[ps.ORGANISM.eq("human") & ~ps.ACC_ID.str.contains("-")].copy()
psh["rsd"] = psh.MOD_RSD.str.extract(r"^([STY]\d+)")[0]
psh = psh[psh.rsd.notna()]
psh["site"] = psh.ACC_ID + "_" + psh.rsd
grp2site = psh.drop_duplicates("SITE_GRP_ID").set_index("SITE_GRP_ID").site
for c in ["LT_LIT", "MS_LIT", "MS_CST"]:
    psh[c] = pd.to_numeric(psh[c], errors="coerce").fillna(0)
siteinfo = psh.groupby("site")[["LT_LIT", "MS_LIT", "MS_CST"]].max()
# map KS rows -> (kinase symbol, human site)
U = pd.read_csv("/workspace/scripts/infra_compendium/kinact_mvp_pkg/model/kinases.tsv", sep="\t", keep_default_na=False)
by_acc = dict(zip(U.uniprot, U.kinase)); by_name = {k.upper(): k for k in U.kinase}
ks = ks[ks.KIN_ORGANISM.isin({"human", "mouse", "rat"}) & ks.SUB_ORGANISM.isin({"human", "mouse", "rat"})].copy()
ks["kinase"] = ks.KIN_ACC_ID.str.split("-").str[0].map(by_acc).where(ks.KIN_ORGANISM.eq("human")).fillna(ks.GENE.str.upper().map(by_name))
rsd = ks.SUB_MOD_RSD.str.extract(r"^([STY]\d+)")[0]
hum = ks.SUB_ORGANISM.eq("human") & ~ks.SUB_ACC_ID.str.contains("-") & rsd.notna()
ks["site"] = (ks.SUB_ACC_ID + "_" + rsd).where(hum).fillna(ks.SITE_GRP_ID.map(grp2site))
ks = ks.dropna(subset=["kinase", "site"])
ks["vivo"] = ks.IN_VIVO_RXN.str.strip().eq("X"); ks["vitro"] = ks.IN_VITRO_RXN.str.strip().eq("X")
ks["hk"] = ks.KIN_ORGANISM.eq("human"); ks["hs"] = ks.SUB_ORGANISM.eq("human")
EA = ks.groupby(["kinase", "site"]).agg(in_vivo=("vivo", "max"), in_vitro=("vitro", "max"), human_kinase_record=("hk", "max"),
                                         human_sub_record=("hs", "max"), n_records=("vivo", "size"),
                                         window=("SITE_+/-7_AA", "first")).reset_index()
E = lib.merge(EA, on=["kinase", "site"], how="left")
E = E.join(siteinfo, on="site")
# regulatory-site references
rg = read_psp(f"{D}/Regulatory_sites")
rg = rg[rg.ORGANISM.eq("human") & rg.MOD_RSD.str.endswith("-p")]
rg["site"] = rg.ACC_ID + "_" + rg.MOD_RSD.str.replace("-p", "", regex=False)
rg = rg.groupby("site").agg(reg_pmids=("PMIDs", lambda x: ";".join(sorted({p.strip() for v in x for p in str(v).split(";") if p.strip()}))),
                            reg_function=("ON_FUNCTION", "first"))
E = E.join(rg, on="site")
E["n_reg_pmids"] = E.reg_pmids.fillna("").apply(lambda s: len([p for p in s.split(";") if p]))
E["evidence_type"] = np.select([E.in_vivo & E.in_vitro, E.in_vivo, E.in_vitro], ["vivo+vitro", "vivo_only", "vitro_only"], "unknown")
SA = load.sites().set_index("site")
E["gene"] = SA.gene.reindex(E.site).values; E["residue"] = E.site.str.split("_").str[1].str[0]
E["in_compendium"] = E.site.isin(SA.index)
E["group"] = E.kinase.map(k2g)
KA = load.kinases().set_index("kinase")
E["manning_group"] = KA.manning_group.reindex(E.kinase).values; E["manning_family"] = KA.manning_family.reindex(E.kinase).values
mot = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
E["motif_pct"] = [mot.at[s, k] if (k in mot.columns and s in mot.index) else np.nan for k, s in zip(E.kinase, E.site)]
E.to_parquet(f"{OUT}/psp_edges_annot.parquet")
print(E.shape, E.group.notna().sum(), E[E.group.notna()].in_compendium.sum())
print(E[E.group.notna()].evidence_type.value_counts().to_dict())

# ---------------- robust z for LoF conditions
conds = sorted(GC.conditionID.unique())
X = load.sites_matrix(conds)
med = X.median(axis=0); mad = (X - med).abs().median(axis=0) * 1.4826
Z = ((X - med) / mad.clip(lower=0.05)).astype(np.float32)
sg = GC.drop_duplicates("conditionID").set_index("conditionID").sign
for c in sg.index[sg < 0]:
    Z[c] = -Z[c]
Z = Z.clip(-8, 8)
Z = Z[Z.notna().any(axis=1)]
Z.reset_index().to_parquet(f"{OUT}/Z_lof.parquet")
print("Z", Z.shape)
