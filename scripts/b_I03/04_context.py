"""I03 step 4: labels, external edges, kinaseSub predictions, motif, PPI, site annotations.
Outputs (cache/b_I03/):
  psp_edges.parquet     kinase, site  (PSP any kinase organism, mapped to human site via SITE_GRP_ID) + flag human_direct
  ext_edges.parquet     kinase, site, src, n_src  (OmniPath phosphorylation edges + iPTMnet; PSP-derived sources dropped,
                        and any pair present in PSP removed) ; tiers: curated (SIGNOR/HPRD/phosphoELM/dbPTM/Reactome/NCI-PID/
                        Li2012/UniProt/PRO) vs textmining (REACH/Sparser/RLIMS-P/BEL)
  known_any.parquet     every kinase-site pair in PSP, OmniPath (all sources) or iPTMnet (for novelty filtering)
  ks_pairs.parquet      kinaseSub release: kinase, site, ks_pct, ks_prob, ks_assign, ks_rank_site
  ppi_pairs.parquet     gene pairs with physical interaction (BioPlex v4 + 293T/HCT116, BioGRID physical, STRING physical exp>=400)
  site_annot.parquet    site, gene, ochoa functional score, PSP regulatory flag/function
"""
import os, sys, zipfile, gzip
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/cache/b_I03"
EXT = "<HACKATHON_ROOT>/kinaseSubMain/external"
REL = "<HACKATHON_ROOT>/kinaseSubMain/results/release"
S = load.sites()
SITES = set(S.site)
K = load.kinases()
kin_genes = set(K.iloc[:, 0]) if "gene" not in K.columns else set(K.gene)

# ---- PSP via site group
ps = pd.read_csv("/workspace/data/Phosphorylation_site_dataset", sep="\t", skiprows=3, low_memory=False, encoding="latin-1",
                 usecols=["ACC_ID", "MOD_RSD", "ORGANISM", "SITE_GRP_ID"])
ps = ps[(ps.ORGANISM == "human") & ~ps.ACC_ID.str.contains("-", na=True)]
ps["site"] = ps.ACC_ID + "_" + ps.MOD_RSD.str.replace("-p", "", regex=False)
g2site = ps.groupby("SITE_GRP_ID").site.apply(set)
ks = pd.read_csv("/workspace/data/Kinase_Substrate_Dataset", sep="\t", skiprows=3, low_memory=False, encoding="latin-1")
rows = []
for r in ks.itertuples():
    kg = str(r.GENE).upper()
    if r.SUB_ORGANISM == "human" and "-" not in str(r.SUB_ACC_ID):
        rows.append((kg, f"{r.SUB_ACC_ID}_{r.SUB_MOD_RSD}", r.KIN_ORGANISM == "human"))
    for s in g2site.get(r.SITE_GRP_ID, ()):
        rows.append((kg, s, False))
PSP = pd.DataFrame(rows, columns=["kinase", "site", "human_direct"]).groupby(["kinase", "site"]).human_direct.max().reset_index()
hm = pd.read_csv("/workspace/results/infra_compendium/psp_ks_edges_human_mapped.tsv", sep="\t")
PSP = pd.concat([PSP, hm.assign(human_direct=True)]).groupby(["kinase", "site"]).human_direct.max().reset_index()
PSP.to_parquet(f"{OUT}/psp_edges.parquet"); print("PSP pairs", len(PSP))
psp_set = set(zip(PSP.kinase, PSP.site))

# ---- External
om = pd.read_csv(f"{EXT}/omnipath/omnipath_enzsub.tsv", sep="\t", low_memory=False)
om = om[om.modification == "phosphorylation"].copy()
om["site"] = om.substrate.astype(str) + "_" + om.residue_type + om.residue_offset.astype(str)
CUR = {"SIGNOR", "SIGNOR_ProtMapper", "HPRD", "HPRD_MIMP", "phosphoELM", "phosphoELM_MIMP", "dbPTM", "Reactome_ProtMapper",
       "NCI-PID_ProtMapper", "Li2012"}
TXT = {"REACH_ProtMapper", "Sparser_ProtMapper", "RLIMS-P_ProtMapper", "BEL-Large-Corpus_ProtMapper"}
erows, allk = [], []
for r in om.itertuples():
    src = set(r.sources.split(";"))
    allk.append((str(r.enzyme_genesymbol).upper(), r.site))
    cur = src & CUR; txt = src & TXT
    if cur:
        erows.append((str(r.enzyme_genesymbol).upper(), r.site, "omnipath:" + "|".join(sorted(cur)), "curated"))
    elif txt:
        erows.append((str(r.enzyme_genesymbol).upper(), r.site, "omnipath:" + "|".join(sorted(txt)), "textmining"))
ip = pd.read_csv(f"{EXT}/iptmnet/ptm.txt", sep="\t", header=None, low_memory=False)
ip = ip[(ip[0] == "PHOSPHORYLATION") & ip[4].astype(str).str.contains("Human") & ip[7].notna()]
for r in ip.itertuples(index=False):
    site = f"{r[2]}_{r[5]}"
    kg = str(r[7]).upper()
    allk.append((kg, site))
    tier = "textmining" if r[1] in ("rlim+", "rlimsp") else "curated"
    erows.append((kg, site, "iptmnet:" + str(r[1]), tier))
EX = pd.DataFrame(erows, columns=["kinase", "site", "src", "tier"])
EX = EX[[(k, s) not in psp_set for k, s in zip(EX.kinase, EX.site)]]
EX = EX.groupby(["kinase", "site"]).agg(src=("src", lambda x: ";".join(sorted(set(x)))),
                                        tier=("tier", lambda x: "curated" if "curated" in set(x) else "textmining")).reset_index()
EX.to_parquet(f"{OUT}/ext_edges.parquet"); print("external non-PSP pairs", len(EX), EX.tier.value_counts().to_dict(),
                                                  "in compendium", EX.site.isin(SITES).sum())
KA = pd.DataFrame(allk, columns=["kinase", "site"])
KA = pd.concat([KA, PSP[["kinase", "site"]]]).drop_duplicates()
KA.to_parquet(f"{OUT}/known_any.parquet")

# ---- kinaseSub
a = pd.read_csv(f"{REL}/site_top5.tsv.gz", sep="\t")
a = a[a.site_id.isin(SITES)]
b = pd.read_csv(f"{REL}/kinase_top200.tsv.gz", sep="\t")
b = b[b.site_id.isin(SITES)]
KS = pd.concat([a[["gene", "site_id", "pct", "prob", "assignment_prob", "rank"]].rename(columns={"rank": "ks_rank_site"}),
                b[["gene", "site_id", "pct", "prob", "assignment_prob"]].assign(ks_rank_site=np.nan)])
KS = KS.groupby(["gene", "site_id"]).agg(ks_pct=("pct", "max"), ks_prob=("prob", "max"), ks_assign=("assignment_prob", "max"),
                                         ks_rank_site=("ks_rank_site", "min")).reset_index().rename(columns={"gene": "kinase", "site_id": "site"})
KS.to_parquet(f"{OUT}/ks_pairs.parquet"); print("kinaseSub pairs", len(KS))

# ---- PPI
pp = []
for f in ["BioPlex_interactionList_v4a.tsv", "BioPlex_293T_Network_10K_Dec_2019.tsv", "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"]:
    d = pd.read_csv(f"{EXT}/bioplex/{f}", sep="\t")
    pp.append(d[["SymbolA", "SymbolB"]].set_axis(["a", "b"], axis=1).assign(src="bioplex"))
z = zipfile.ZipFile(f"{EXT}/biogrid/BIOGRID-ORGANISM-LATEST.tab3.zip")
nm = [n for n in z.namelist() if "Homo_sapiens" in n][0]
bg = pd.read_csv(z.open(nm), sep="\t", usecols=["Official Symbol Interactor A", "Official Symbol Interactor B",
                                                   "Experimental System Type", "Organism ID Interactor A", "Organism ID Interactor B"], low_memory=False)
bg = bg[(bg["Experimental System Type"] == "physical") & (bg["Organism ID Interactor A"] == 9606) & (bg["Organism ID Interactor B"] == 9606)]
pp.append(bg.iloc[:, :2].set_axis(["a", "b"], axis=1).assign(src="biogrid"))
info = pd.read_csv(f"{EXT}/string/9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1]); info.columns = ["e", "g"]
e2g = dict(zip(info.e, info.g))
st = pd.read_csv(f"{EXT}/string/9606.protein.physical.links.detailed.v12.0.txt.gz", sep=" ", usecols=["protein1", "protein2", "experimental"])
st = st[st.experimental >= 400]
pp.append(pd.DataFrame(dict(a=st.protein1.map(e2g), b=st.protein2.map(e2g), src="string_phys")))
PP = pd.concat(pp).dropna()
PP["a"] = PP.a.astype(str).str.upper(); PP["b"] = PP.b.astype(str).str.upper()
PP2 = pd.concat([PP, PP.rename(columns={"a": "b", "b": "a"})])
PP2 = PP2.groupby(["a", "b"]).src.apply(lambda x: ";".join(sorted(set(x)))).reset_index()
PP2.to_parquet(f"{OUT}/ppi_pairs.parquet"); print("PPI pairs", len(PP2))

# ---- site annotation
fs = pd.read_excel(f"{EXT}/ochoa2020/41587_2019_344_MOESM5_ESM.xlsx", sheet_name="functional_score")
res = S.set_index(["uniprot", "position"]).residue
fs["site"] = fs.uniprot + "_" + [res.get((u, p), "X") for u, p in zip(fs.uniprot, fs.position)] + fs.position.astype(str)
A = S[["site", "uniprot", "residue", "position", "gene", "n_conditions", "n_studies", "psp_regulatory", "psp_reg_function", "sequence_window"]]
A = A.merge(fs[["site", "functional_score"]].drop_duplicates("site"), on="site", how="left")
A.to_parquet(f"{OUT}/site_annot.parquet"); print("functional score coverage", A.functional_score.notna().mean())
