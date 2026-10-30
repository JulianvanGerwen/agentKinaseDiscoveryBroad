"""Extract kinaseSubMain release percentile scores (kinase x site) for compendium sites measured in >=20 conditions.
Output: cache/pct.parquet (sites x kinases, float32 percentile within kinase over U_obs; NaN = residue mismatch)."""
import numpy as np, pandas as pd, pyarrow.parquet as pq
from pathlib import Path
R = "<HACKATHON_ROOT>/kinaseSubMain"
OUT = Path("/workspace/results/b_I01/cache"); OUT.mkdir(parents=True, exist_ok=True)
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", usecols=["site", "n_conditions", "n_studies"])
S = S[S.n_conditions >= 20]
ks = pd.read_parquet(f"{R}/cache/index/kinases.parquet")
si = pd.read_parquet(f"{R}/cache/index/sites.parquet")[["site_id"]].reset_index(drop=True)
si["site_idx"] = np.arange(len(si))
si = si[si.site_id.isin(S.site)]
keep = np.zeros(350736, bool); keep[si.site_idx.values] = True
pos = -np.ones(350736, np.int64); pos[si.site_idx.values] = np.arange(len(si))
M = np.full((len(si), len(ks)), np.nan, np.float32)
f = pq.ParquetFile(f"{R}/cache/release/r2r3_release__20261027-012616__s0/release_scores.parquet")
for g in range(f.num_row_groups):
    t = f.read_row_group(g, columns=["kinase_idx", "site_idx", "pct"]).to_pandas()
    t = t[keep[t.site_idx.values]]
    M[pos[t.site_idx.values], t.kinase_idx.values.astype(int)] = t.pct.values
df = pd.DataFrame(M, index=si.site_id.values, columns=ks.gene.values if "gene" in ks else ks.iloc[:, 1].values)
df.to_parquet(OUT / "pct.parquet")
print(df.shape, np.isnan(M).mean())
