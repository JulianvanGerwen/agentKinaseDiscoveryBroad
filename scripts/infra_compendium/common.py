"""Shared paths / helpers for the infra_compendium build (all outputs under /workspace)."""
import os
from io import StringIO
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import pandas as pd

WS = Path("/workspace")
DATA = WS / "data"
RES = WS / "results" / "infra_compendium"
CACHE = WS / "cache" / "compendium"
SCRIPTS = WS / "scripts" / "infra_compendium"
KAM = Path("<HACKATHON_ROOT>/kinaseActivityMain")
KSM = Path("<HACKATHON_ROOT>/kinaseSubMain")
for p in (RES, CACHE):
    p.mkdir(parents=True, exist_ok=True)


def read_psp(path) -> pd.DataFrame:
    """PSP flat files: header lines before the 'GENE' / column line, latin-1, TSV."""
    with open(path, "rt", encoding="latin-1") as fh:
        lines = fh.readlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(("GENE\t", "KINASE\t")))
    return pd.read_csv(StringIO("".join(lines[start:])), sep="\t", dtype=str, keep_default_na=False)


def kinase_universe() -> pd.DataFrame:
    """503-kinase universe of core-v3 (gene symbol, uniprot, Manning group/family)."""
    return pd.read_csv(SCRIPTS / "kinact_mvp_pkg" / "model" / "kinases.tsv", sep="\t", keep_default_na=False)
