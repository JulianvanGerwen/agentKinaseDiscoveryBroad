"""Score all compendium conditions with core-v3 (MVP copy + evaluated-equivalent edges).

Usage: python 03_score_corev3.py [n_conditions_for_timing_test]
Output: cache/compendium/corev3_long.parquet (condition, kinase, z, z_core, z_sig, mean, sd, n_substrates)
"""
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "8")
import pandas as pd
import pyarrow.parquet as pq
from threadpoolctl import threadpool_limits

from common import CACHE, RES, SCRIPTS

sys.path.insert(0, str(SCRIPTS / "kinact_mvp_pkg"))
from kinact_mvp import KinaseActivityModel  # noqa: E402

N = int(sys.argv[1]) if len(sys.argv) > 1 else None
t0 = time.time()
model = KinaseActivityModel.from_config(SCRIPTS / "kinact_mvp_pkg/model/config.json")
print(f"model loaded {time.time()-t0:.1f}s", flush=True)
cols = pq.read_schema(RES / "site_matrix.parquet").names[1:]
cols = cols[:N] if N else cols
out = []
CH = 250
for i in range(0, len(cols), CH):
    cc = cols[i:i + CH]
    X = pq.read_table(RES / "site_matrix.parquet", columns=["site"] + cc).to_pandas().set_index("site")
    t1 = time.time()
    with threadpool_limits(limits=8):
        res = model.score(X, device="cpu", n_jobs=8)
    lg = res.long()
    out.append(lg)
    print(f"{i+len(cc)}/{len(cols)} conditions, {time.time()-t1:.1f}s (total {time.time()-t0:.0f}s)", flush=True)
df = pd.concat(out, ignore_index=True)
print(df.head(), df.columns.tolist())
suffix = f"_test{N}" if N else ""
df.to_parquet(CACHE / f"corev3_long{suffix}.parquet", index=False)
print(f"done {time.time()-t0:.0f}s; info {res.info}")
