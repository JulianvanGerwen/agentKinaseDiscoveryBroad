"""Orthogonal kinase x kinase label matrices (503 kinases): DepMap co-dependency, GTEx co-expression,
PPI (experimental screens: BioPlex/OpenCell/HuRI/HumanCellMap; literature: BioGRID physical/STRING physical>=700),
co-localisation (HPA Jaccard; COMPARTMENTS cosine)."""
import numpy as np, pandas as pd, zipfile, gzip, io, sys
from pathlib import Path
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
E = Path("<HACKATHON_ROOT>/kinaseSubMain/external")
CA = Path("/workspace/results/b_I01/cache")
K = load.kinases(); kin = K.kinase.tolist(); ks = set(kin)
H = pd.read_csv(E / "hgnc/hgnc_complete_set.txt", sep="\t", low_memory=False, usecols=["symbol", "prev_symbol", "alias_symbol", "ensembl_gene_id", "entrez_id", "uniprot_ids"])
sym = {}
for _, r in H.iterrows():
    if r.symbol in ks:
        for c in ["prev_symbol", "alias_symbol"]:
            if isinstance(r[c], str):
                for a in r[c].split("|"): sym.setdefault(a, r.symbol)
for k in kin: sym[k] = k
ens = {r.ensembl_gene_id: r.symbol for _, r in H.iterrows() if isinstance(r.ensembl_gene_id, str) and r.symbol in ks}
up = {}
for _, r in K.iterrows(): up[r.uniprot] = r.kinase
m = lambda g: sym.get(g)
out = {}

def pairs_to_mat(pairs, name):
    A = pd.DataFrame(0, index=kin, columns=kin, dtype=np.int8)
    n = 0
    for a, b in pairs:
        if a and b and a != b and a in ks and b in ks:
            A.loc[a, b] = 1; A.loc[b, a] = 1; n += 1
    print(name, "kinase-kinase edges", int(A.values.sum() // 2))
    return A

# DepMap
D = pd.read_csv(E / "depmap/DepMap_26Q1_chronos_gene_effect.csv", index_col=0)
D.columns = [c.split(" (")[0] for c in D.columns]
D = D.loc[:, [c for c in D.columns if m(c)]]; D.columns = [m(c) for c in D.columns]
D = D.loc[:, ~D.columns.duplicated()]
print("depmap kinases", D.shape)
dep = D.corr(min_periods=200).reindex(index=kin, columns=kin)
dep.to_parquet(CA / "lab_depmap.parquet")
# GTEx
G = pd.read_csv(E / "gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
G["sym"] = G.Description.map(m); G = G.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym").iloc[:, 2:]
G = np.log2(G.astype(float) + 1); G = G[(G.max(axis=1) > 1)]
gt = G.T.rank().corr().reindex(index=kin, columns=kin)
gt.to_parquet(CA / "lab_gtex.parquet")
# PPI experimental
pe = []
for f in ["BioPlex_293T_Network_10K_Dec_2019.tsv", "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv", "BioPlex_interactionList_v4a.tsv"]:
    b = pd.read_csv(E / "bioplex" / f, sep="\t")
    pe += [(up.get(a.split("-")[0]) or m(sa), up.get(c.split("-")[0]) or m(sb)) for a, c, sa, sb in zip(b.UniprotA, b.UniprotB, b.SymbolA, b.SymbolB)]
oc = pd.read_csv(E / "opencell/opencell-protein-interactions.csv")
pe += [(m(a), m(b)) for a, b in zip(oc.target_gene_name, oc.interactor_gene_name)]
hu = pd.read_csv(E / "huri/HI-union.tsv", sep="\t", header=None)
pe += [(ens.get(a), ens.get(b)) for a, b in zip(hu[0], hu[1])]
hc = pd.read_csv(E / "humancellmap/saint-latest.txt", sep="\t", usecols=["Bait", "PreyGene", "BFDR"])
hc = hc[hc.BFDR <= 0.01]
pe += [(m(a), m(b)) for a, b in zip(hc.Bait, hc.PreyGene)]
ppi_exp = pairs_to_mat(pe, "ppi_exp")
# PPI literature
z = zipfile.ZipFile(E / "biogrid/BIOGRID-ORGANISM-LATEST.tab3.zip")
nm = [x for x in z.namelist() if "Homo_sapiens" in x][0]
bg = pd.read_csv(z.open(nm), sep="\t", usecols=["Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System Type", "Organism ID Interactor A", "Organism ID Interactor B"], low_memory=False)
bg = bg[(bg["Experimental System Type"] == "physical") & (bg["Organism ID Interactor A"] == 9606) & (bg["Organism ID Interactor B"] == 9606)]
pl = [(m(a), m(b)) for a, b in zip(bg["Official Symbol Interactor A"], bg["Official Symbol Interactor B"])]
info = pd.read_csv(E / "string/9606.protein.info.v12.0.txt.gz", sep="\t")
smap = dict(zip(info.iloc[:, 0], info.preferred_name.map(m)))
st = pd.read_csv(E / "string/9606.protein.physical.links.detailed.v12.0.txt.gz", sep=" ")
st = st[st.combined_score >= 700]
pl += [(smap.get(a), smap.get(b)) for a, b in zip(st.protein1, st.protein2)]
ppi_lit = pairs_to_mat(pl, "ppi_lit")
ppi_exp.to_parquet(CA / "lab_ppi_exp.parquet"); ppi_lit.to_parquet(CA / "lab_ppi_lit.parquet")
# HPA co-localisation (Jaccard of main+additional locations)
z = zipfile.ZipFile(E / "hpa/subcellular_location.tsv.zip")
hp = pd.read_csv(z.open("subcellular_location.tsv"), sep="\t")
hp["sym"] = hp["Gene name"].map(m); hp = hp.dropna(subset=["sym"]).drop_duplicates("sym")
locs = {r.sym: set(str(r["Main location"]).split(";")) | (set(str(r["Additional location"]).split(";")) if isinstance(r["Additional location"], str) else set()) for _, r in hp.iterrows()}
J = pd.DataFrame(np.nan, index=kin, columns=kin)
for a in locs:
    for b in locs:
        J.loc[a, b] = len(locs[a] & locs[b]) / len(locs[a] | locs[b])
J.to_parquet(CA / "lab_hpa.parquet")
# COMPARTMENTS integrated (score>=2) cosine
cp = pd.read_csv(E / "compartments/human_compartment_integrated_full.tsv", sep="\t", header=None, names=["id", "name", "go", "term", "score"])
cp = cp[cp.score >= 2]; cp["sym"] = cp.name.map(m); cp = cp.dropna(subset=["sym"])
V = cp.pivot_table(index="sym", columns="go", values="score", aggfunc="max").fillna(0)
V = V.loc[:, (V > 0).sum() >= 5]
Vn = V.div(np.sqrt((V ** 2).sum(axis=1)), axis=0)
cos = (Vn @ Vn.T).reindex(index=kin, columns=kin)
cos.to_parquet(CA / "lab_compartments.parquet")
for n, A in [("depmap", dep), ("gtex", gt), ("hpa", J), ("comp", cos)]:
    print(n, "kinases covered", int(A.notna().any(axis=1).sum()))
