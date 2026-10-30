"""Site-level sequence features for phosphatase preference (compendium + CPTAC sites)."""
import gzip, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
OUT = "/workspace/results/b_I13/cache"
FA = "<HACKATHON_ROOT>/kinaseSubMain/external/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz"
seqs, cur, buf = {}, None, []
for line in gzip.open(FA, "rt"):
    if line.startswith(">"):
        if cur: seqs[cur] = "".join(buf)
        cur = line.split("|")[1]; buf = []
    else: buf.append(line.strip())
seqs[cur] = "".join(buf)

S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", usecols=["site", "uniprot", "residue", "position", "gene"])
keys = set(S.site)
for c in cl.COHORTS:
    keys |= set(pd.read_parquet(f"/workspace/results/infra_cptac/cohorts/{c}/phospho.parquet", columns=[]).index)
keys = sorted(keys)
df = pd.DataFrame({"site": keys})
sp = df.site.str.rsplit("_", n=1, expand=True)
df["uniprot"] = sp[0]; df["residue"] = sp[1].str[0]; df["position"] = pd.to_numeric(sp[1].str[1:], errors="coerce")
df = df[df.uniprot.isin(seqs) & df.position.notna()].copy(); df.position = df.position.astype(int)

# protein-level SLiMs
SLIM = {"LxxIxE": r"(?=([LMFI]..[IVL].E))",          # PP2A-B56 (Hertz 2016)
        "RVxF": r"(?=([RK].{0,1}[VI][^P][FW]))",      # PP1 docking (consensus, permissive)
        "PxIxIT": r"(?=(P.I.[IV][TSHDE]))",          # calcineurin
        "LxVP": r"(?=(L.VP))"}                        # calcineurin
slim_pos = {k: {} for k in SLIM}
for u in df.uniprot.unique():
    s = seqs[u]
    for k, p in SLIM.items(): slim_pos[k][u] = np.array([m.start() + 1 for m in re.finditer(p, s)])

def win(u, p, f=7):
    s = seqs[u]; 
    return "".join(s[i] if 0 <= i < len(s) else "_" for i in range(p - 1 - f, p + f))
df["window"] = [win(u, p) for u, p in zip(df.uniprot, df.position)]
df = df[df.window.str[7] == df.residue].copy()   # residue check
W = np.array([list(w) for w in df.window])        # n x 15, centre=7
def cnt(cols, aas): return np.isin(W[:, cols], list(aas)).sum(1)
L = lambda a, b: [7 + i for i in range(a, b + 1)]
df["isT"] = (df.residue == "T").astype(int); df["isY"] = (df.residue == "Y").astype(int)
df["P1"] = (W[:, 8] == "P").astype(int)
df["Q1"] = (W[:, 8] == "Q").astype(int)
df["R_m3"] = np.isin(W[:, 4], ["R", "K"]).astype(int)
df["R_m2"] = np.isin(W[:, 5], ["R", "K"]).astype(int)
df["basic_m7m1"] = cnt(L(-7, -1), "KR"); df["basic_p1p7"] = cnt(L(1, 7), "KR")
df["basic_m3p3"] = cnt(L(-3, -1) + L(1, 3), "KR")
df["acid_m3p3"] = cnt(L(-3, -1) + L(1, 3), "DE")
df["acid_p1p3"] = cnt(L(1, 3), "DE")
df["acid_m7p7"] = cnt(L(-7, -1) + L(1, 7), "DE")
df["pro_m7p7"] = cnt(L(-7, -1) + L(2, 7), "P")
df["hydrophob_p1"] = np.isin(W[:, 8], list("LIVMF")).astype(int)
df["ST_m7p7"] = cnt(L(-7, -1) + L(1, 7), "ST")
df["bipartite_basic"] = ((df.basic_m7m1 > 0) & (df.basic_p1p7 > 0)).astype(int)
for k in SLIM:
    has, dist = [], []
    for u, p in zip(df.uniprot, df.position):
        a = slim_pos[k][u]; has.append(len(a)); dist.append(np.abs(a - p).min() if len(a) else 5000)
    df[f"n_{k}"] = has; df[f"logdist_{k}"] = np.log10(np.array(dist) + 1)
    df[f"has_{k}"] = (np.array(has) > 0).astype(int)
df["prot_len_log"] = np.log10([len(seqs[u]) for u in df.uniprot])
# literature rule B55 score (Cundell 2016 / Kruse 2020 / Hein 2017): pT, S/T-P, basic residues both sides
df["B55_lit"] = 1.0 * df.isT + 1.0 * df.P1 + 0.5 * np.minimum(df.basic_m7m1, 3) + 0.5 * np.minimum(df.basic_p1p7, 3) - 0.5 * np.minimum(df.acid_m3p3, 3)
df["B56_lit"] = 1.0 * df.R_m3 + 1.0 * (df.logdist_LxxIxE < 2.0) + 0.5 * df.isT * 0  # basophilic + nearby LxxIxE (<100 aa)
df = df[df.residue.isin(list("STY"))]
df.to_parquet(f"{OUT}/site_features.parquet")
print(df.shape, df.residue.value_counts().to_dict()); print(df.describe().T.iloc[:, :3])
