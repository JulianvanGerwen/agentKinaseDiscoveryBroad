"""Site/protein annotation for b_I17: PSP + OmniPath kinase-substrate edges, kinaseSub predictions (top5),
STRING (combined>=700) graph, COMPARTMENTS main localisation, DepMap mean Chronos.
Outputs cache/b_I17/{ks_edges.parquet, pred_top5.parquet, string_edges.parquet, uniprot_gene.tsv, loc.tsv, depmap_mean.tsv}
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import gzip, numpy as np, pandas as pd, sys
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
EXT = "<HACKATHON_ROOT>/kinaseSubMain/external"
REL = "<HACKATHON_ROOT>/kinaseSubMain/results/release"
OUT = "/workspace/cache/b_I17"
S = load.sites()[["site", "uniprot", "residue", "position", "gene", "psp_kinases_human", "pred_top1_kinase", "pred_top3_kinases"]]
S.to_parquet(f"{OUT}/sites_min.parquet")
u2g = S.dropna(subset=["gene"]).drop_duplicates("uniprot").set_index("uniprot").gene

# kinase-substrate edges (site level): PSP (from sites.tsv) + OmniPath enzsub
ks = []
for s, k in S.dropna(subset=["psp_kinases_human"])[["site", "psp_kinases_human"]].itertuples(index=False):
    for kk in str(k).split(";"):
        if kk: ks.append((s, kk.strip(), "PSP"))
om = pd.read_csv(f"{EXT}/omnipath/omnipath_enzsub.tsv", sep="\t")
om = om[om.modification == "phosphorylation"]
om["site"] = om.substrate + "_" + om.residue_type + om.residue_offset.astype(str)
ks += list(zip(om.site, om.enzyme_genesymbol, ["OmniPath"] * len(om)))
KS = pd.DataFrame(ks, columns=["site", "kinase", "src"]).drop_duplicates()
KS.to_parquet(f"{OUT}/ks_edges.parquet")
# protein-level directed enzyme->substrate (for network distance)
pe = om[["enzyme_genesymbol", "substrate_genesymbol"]].drop_duplicates()
pe.columns = ["a", "b"]; pe.to_parquet(f"{OUT}/omnipath_protein_edges.parquet")

P = pd.read_csv(f"{REL}/site_top5.tsv.gz", sep="\t")
P = P[P.site_id.isin(set(S.site))]
P.to_parquet(f"{OUT}/pred_top5.parquet")

# STRING
al = pd.read_csv(f"{EXT}/string/9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1])
al.columns = ["ensp", "gene"]
e2g = dict(zip(al.ensp, al.gene))
L = pd.read_csv(f"{EXT}/string/9606.protein.links.detailed.v12.0.txt.gz", sep=" ", usecols=["protein1", "protein2", "experimental", "database", "combined_score"])
L = L[L.combined_score >= 700]
L["a"] = L.protein1.map(e2g); L["b"] = L.protein2.map(e2g)
L[["a", "b", "combined_score"]].to_parquet(f"{OUT}/string_edges.parquet")

# COMPARTMENTS integrated: main locs (score>=3)
cp = pd.read_csv(f"{EXT}/compartments/human_compartment_integrated_full.tsv", sep="\t", header=None, names=["id", "gene", "go", "term", "score"])
locs = {"Nucleus": "GO:0005634", "Cytosol": "GO:0005829", "Plasma membrane": "GO:0005886", "Mitochondrion": "GO:0005739"}
cp = cp[cp.go.isin(locs.values()) & (cp.score >= 3)]
loc = cp.pivot_table(index="gene", columns="term", values="score", aggfunc="max")
loc.to_csv(f"{OUT}/loc.tsv", sep="\t")

# DepMap mean chronos (no Model.csv available -> pan-line mean and fraction < -0.5)
D = pd.read_csv(f"{EXT}/depmap/DepMap_26Q1_chronos_gene_effect.csv", index_col=0)
D.columns = [c.split(" (")[0] for c in D.columns]
dm = pd.DataFrame(dict(mean_chronos=D.mean(), frac_dep=(D < -0.5).mean()))
dm.to_csv(f"{OUT}/depmap_mean.tsv", sep="\t")
print(len(KS), len(P), len(L), loc.shape, dm.shape)
