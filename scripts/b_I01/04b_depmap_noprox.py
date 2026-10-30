"""DepMap co-dependency with chromosomal-proximity pairs removed (same chromosome arm + major band, e.g. 19q13),
because CRISPR (Chronos) gene effects of neighbouring genes correlate through copy-number/cutting artefacts."""
import pandas as pd, numpy as np, re
from pathlib import Path
E = Path("<HACKATHON_ROOT>/kinaseSubMain/external"); CA = Path("/workspace/results/b_I01/cache")
D = pd.read_parquet(CA / "lab_depmap.parquet")
H = pd.read_csv(E / "hgnc/hgnc_complete_set.txt", sep="\t", low_memory=False, usecols=["symbol", "location"]).set_index("symbol")
def band(x):
    m = re.match(r"^(\d+|X|Y)([pq]\d)", str(x)); return m.group(0) if m else None
b = pd.Series({k: band(H.location.get(k)) for k in D.index}, dtype=object)
bv = np.array(b.tolist(), dtype=object)
same = (bv[:, None] == bv[None, :]) & b.notna().values[:, None]
print("pairs removed", same.sum() // 2, "kinases with band", b.notna().sum())
D.where(~same).to_parquet(CA / "lab_depmap_noprox.parquet")
b.to_csv(CA / "kinase_band.tsv", sep="\t")
