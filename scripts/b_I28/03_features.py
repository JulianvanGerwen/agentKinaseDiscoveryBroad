"""I28 step 3: site/protein features for all atlas sites (substrates + random controls).
Out: cache/b_I28/features.parquet (index site)
"""
import os, sys, re, gzip
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from lib28 import *
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
EXT = f"{KSM}/external"
L = pd.read_parquet(f"{C28}/L_units.parquet")
sites = sorted(L.site.unique())
SA = load.sites().set_index("site")
F = pd.DataFrame(index=pd.Index(sites, name="site"))
F["acc"] = F.index.str.split("_").str[0]; F["pos"] = F.index.str.split("_").str[1].str[1:].astype(int)
F["res"] = F.index.str.split("_").str[1].str[0]
w = SA.sequence_window.reindex(sites).fillna("")
def wpos(s, i):
    return s[7 + i] if len(s) == 15 else ""
F["is_T"] = (F.res == "T").astype(float); F["is_Y"] = (F.res == "Y").astype(float)
F["P_plus1"] = [float(wpos(s, 1) == "P") for s in w]
F["n_basic_pm5"] = [sum(c in "KR" for c in (s[2:7] + s[8:13])) for s in w]
F["n_acidic_pm5"] = [sum(c in "DE" for c in (s[2:7] + s[8:13])) for s in w]
F["basic_minus3"] = [float(wpos(s, -3) in "KR") if s else np.nan for s in w]
F["psp_regulatory"] = SA.psp_regulatory.reindex(sites).astype(float).values
F["log_n_conditions"] = np.log10(SA.n_conditions.reindex(sites).astype(float).values)
F["log_n_studies"] = np.log10(SA.n_studies.reindex(sites).astype(float).values)
# multisite cluster: other compendium phosphosites on same protein within +-5 / +-15
pos_by = SA.groupby("uniprot").position.apply(lambda x: np.sort(x.astype(int).values)).to_dict()
F["n_sites_pm5"] = [((np.abs(pos_by.get(a, np.array([])) - p) <= 5).sum() - 1) for a, p in zip(F.acc, F.pos)]
F["n_sites_pm15"] = [((np.abs(pos_by.get(a, np.array([])) - p) <= 15).sum() - 1) for a, p in zip(F.acc, F.pos)]
# structure
st = pd.read_parquet(f"{KSM}/cache/structure/site_structure.parquet", columns=["site_id", "plddt", "plddt_s10", "rsa", "idr", "high_acc", "ss3"],
                     filters=[("site_id", "in", sites)]).set_index("site_id")
F = F.join(st)
F["helix"] = (F.ss3 == "H").astype(float).where(F.ss3.notna()); F["strand"] = (F.ss3 == "E").astype(float).where(F.ss3.notna())
F = F.drop(columns="ss3"); F["idr"] = F.idr.astype(float); F["high_acc"] = F.high_acc.astype(float)
# Ochoa 2020
oc_c = f"{C28}/ochoa.parquet"
if not os.path.exists(oc_c):
    o = pd.read_excel(f"{EXT}/ochoa2020/41587_2019_344_MOESM4_ESM.xlsx", sheet_name="annotated_phosphoproteome",
                      usecols=["uniprot", "position", "residue", "disopred_score", "w0_mya", "EV_max_frequency", "adj_ptms_w21",
                               "paxdb_abundance_log10", "isInterface", "pubmed_counts", "quant_top1", "isHotspot"])
    fs = pd.read_excel(f"{EXT}/ochoa2020/41587_2019_344_MOESM5_ESM.xlsx")
    o = o.merge(fs, on=["uniprot", "position"], how="left")
    o["site"] = o.uniprot + "_" + o.residue + o.position.astype(str)
    o.drop(columns=["uniprot", "position", "residue"]).drop_duplicates("site").to_parquet(oc_c)
o = pd.read_parquet(oc_c).set_index("site")
for c in ["isInterface", "isHotspot"]:
    o[c] = o[c].astype(float)
F = F.join(o)
F["log_pubmed"] = np.log10(1 + F.pop("pubmed_counts"))
# localisation (COMPARTMENTS integrated, confidence >=3), by gene symbol
comp = pd.read_csv(f"{EXT}/compartments/human_compartment_integrated_full.tsv", sep="\t", header=None, names=["id", "name", "go", "term", "score"])
GO = {"loc_nucleus": "GO:0005634", "loc_cytosol": "GO:0005829", "loc_plasma_membrane": "GO:0005886", "loc_mitochondrion": "GO:0005739",
      "loc_ER": "GO:0005783", "loc_cytoskeleton": "GO:0005856", "loc_nucleolus": "GO:0005730", "loc_chromosome": "GO:0005694"}
gene = SA.gene.reindex(sites)
for k, g in GO.items():
    s = set(comp[(comp.go == g) & (comp.score >= 3)].name)
    F[k] = [float(x in s) if isinstance(x, str) else np.nan for x in gene]
# complex/PPI degree (BioPlex 293T + OpenCell)
bp = pd.read_csv(f"{EXT}/bioplex/BioPlex_293T_Network_10K_Dec_2019.tsv", sep="\t")
deg = pd.concat([bp.UniprotA.str.split("-").str[0], bp.UniprotB.str.split("-").str[0]]).value_counts()
F["log_bioplex_degree"] = np.log10(1 + F.acc.map(deg).fillna(0))
ocl = pd.read_csv(f"{EXT}/opencell/opencell-protein-interactions.csv")
gcol = [c for c in ocl.columns if "gene_name" in c]
odeg = pd.concat([ocl[c] for c in gcol]).value_counts() if gcol else pd.Series(dtype=float)
F["log_opencell_degree"] = np.log10(1 + gene.map(odeg).fillna(0).values)
# phosphatase docking motifs (sequence)
seqs = {}
with gzip.open(f"{EXT}/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz", "rt") as fh:
    acc, buf = None, []
    for line in fh:
        if line.startswith(">"):
            if acc: seqs[acc] = "".join(buf)
            acc = line.split("|")[1]; buf = []
        else:
            buf.append(line.strip())
    if acc: seqs[acc] = "".join(buf)
MOT = {"PP1_RVxF": r"(?=([RK].?[VI][^P][FW]))", "B56_LxxIxE": r"(?=(L..[ILV].E))", "CN_PxIxIT": r"(?=(P.I.[IV][TSHN]))", "CN_LxVP": r"(?=(L.VP))"}
for m, rx in MOT.items():
    dist, dens = [], []
    for a, p in zip(F.acc, F.pos):
        s = seqs.get(a)
        if s is None:
            dist.append(np.nan); dens.append(np.nan); continue
        hits = np.array([x.start() + 1 for x in re.finditer(rx, s)])
        dist.append(np.log10(1 + (np.abs(hits - p).min() if len(hits) else len(s))))
        dens.append(1000 * len(hits) / len(s))
    F[f"logdist_{m}"] = dist; F[f"dens_{m}"] = dens
F["log_prot_len"] = np.log10([len(seqs.get(a, "")) or np.nan for a in F.acc])
F.drop(columns=["acc", "pos", "res"]).to_parquet(f"{C28}/features.parquet")
print(F.shape); print(F.notna().mean().round(2).to_string())
