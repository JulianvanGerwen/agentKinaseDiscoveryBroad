"""Kinase-substrate edge library for core-v3 (written into the /workspace copy of the MVP package).

PRIMARY (used for scoring) = the library the evaluated core-v3 used:
  cur  = kinaseActivityMain cache/libraries/curated_union.parquet (PSP human/mouse/rat + OmniPath curated + iPTMnet,
         with the evaluation leakage filter; 27,598 edges, sha256 in results/libraries/build_manifest.json)
  r2np = mvp/model/predicted_edges.tsv.gz (kinaseSubMain R2-no-phospho top-100 per kinase), minus curated pairs
-> scripts/infra_compendium/kinact_mvp_pkg/model/edges.tsv.gz  (+ copy in cache/compendium/edges_core_v3_eval.tsv.gz)

ALTERNATIVE (documented, not used for scoring): mvp/build_edges.py recipe on data/ PSP 2025-05 + OmniPath 2026-10-26
-> cache/compendium/edges_mvp_rebuild.tsv.gz ; overlap reported.
"""
import hashlib
import shutil
import subprocess
import sys

import pandas as pd

from common import CACHE, DATA, KAM, RES, SCRIPTS

PKG = SCRIPTS / "kinact_mvp_pkg"
MODEL = PKG / "model"

cu = pd.read_parquet(KAM / "cache/libraries/curated_union.parquet")
kin = pd.read_csv(MODEL / "kinases.tsv", sep="\t", keep_default_na=False)
cu = cu[cu["kinase"].isin(set(kin["kinase"]))]
cur = cu[["site", "kinase"]].drop_duplicates().assign(weight=1.0, source="cur")
pred = pd.read_csv(MODEL / "predicted_edges.tsv.gz", sep="\t")
key = set(zip(cur["site"], cur["kinase"]))
pred = pred[[(s, k) not in key for s, k in zip(pred["site"], pred["kinase"])]]
edges = pd.concat([cur, pred[["site", "kinase", "weight", "source"]]], ignore_index=True)
edges = edges.sort_values(["kinase", "source", "site"]).reset_index(drop=True)
out = MODEL / "edges.tsv.gz"
edges.to_csv(out, sep="\t", index=False, compression="gzip")
shutil.copy(out, CACHE / "edges_core_v3_eval.tsv.gz")
sha = hashlib.sha256(out.read_bytes()).hexdigest()
print(f"primary: {len(cur)} curated edges ({cur.kinase.nunique()} kinases), {len(pred)} r2np edges; sha256 {sha[:16]}")

# alternative MVP recipe (PSP from data/, OmniPath copy in kinaseActivityMain/external)
alt = CACHE / "edges_mvp_rebuild.tsv.gz"
subprocess.run([sys.executable, str(PKG / "build_edges.py"), "--psp-ks", str(DATA / "Kinase_Substrate_Dataset"),
                "--psp-sites", str(DATA / "Phosphorylation_site_dataset"),
                "--omnipath", str(KAM / "external/omnipath/omnipath_enzsub.tsv"), "--out", str(alt)], check=True)
a = pd.read_csv(alt, sep="\t")
ac = set(zip(a.loc[a.source == "cur", "site"], a.loc[a.source == "cur", "kinase"]))
print(f"alt curated {len(ac)}; primary curated {len(key)}; shared {len(ac & key)} "
      f"({len(ac & key)/len(key):.3f} of primary, {len(ac & key)/len(ac):.3f} of alt)")
pre = pd.read_csv(MODEL / "edges_repo_prebuilt.tsv.gz", sep="\t")
pc = set(zip(pre.loc[pre.source == "cur", "site"], pre.loc[pre.source == "cur", "kinase"]))
print(f"repo prebuilt mvp/model/edges.tsv.gz curated {len(pc)}; identical to alt rebuild: {pc == ac}")
with open(RES / "edges_provenance.txt", "w") as fh:
    fh.write(f"primary edges.tsv.gz sha256 {sha}\ncurated(curated_union) {len(cur)} edges / {cur.kinase.nunique()} kinases; "
             f"r2np {len(pred)} edges / {pred.kinase.nunique()} kinases\n"
             f"alt MVP rebuild curated {len(ac)}; shared with primary {len(ac & key)}\n"
             f"repo mvp/model/edges.tsv.gz (prebuilt, 2026-10-28) curated {len(pc)}; identical to alt rebuild: {pc == ac}\n")
