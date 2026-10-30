"""Loaders for the infra_compendium assets.  `import sys; sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load`

    X    = load.sites_matrix(conditions=None)       # DataFrame sites x conditions (float32, NaN = not measured)
    S    = load.sites()                             # site annotation (gene, residue, PSP flags, predicted kinases)
    C    = load.conditions()                        # condition annotation
    A    = load.activity("corev3")                  # kinase x condition z (wide); load.activity("psp_z")
    N    = load.activity("corev3", "nsub")          # n measured substrates (wide)
    L    = load.activity_long()                     # long, both methods
    K    = load.kinases(); D = load.drug_targets()
"""
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

RES = Path("/workspace/results/infra_compendium")


def sites_matrix(conditions=None, sites=None) -> pd.DataFrame:
    cols = None if conditions is None else ["site"] + list(conditions)
    df = pq.read_table(RES / "site_matrix.parquet", columns=cols).to_pandas().set_index("site")
    return df if sites is None else df.loc[df.index.intersection(sites)]


def sites() -> pd.DataFrame:
    return pd.read_csv(RES / "sites.tsv", sep="\t", low_memory=False)


def conditions() -> pd.DataFrame:
    return pd.read_csv(RES / "conditions_annotated.tsv", sep="\t", dtype={"study_pmid": str}, keep_default_na=True)


def activity(method: str = "corev3", what: str = "z") -> pd.DataFrame:
    tag = {"corev3": "corev3", "psp_z": "pspz", "pspz": "pspz"}[method]
    return pd.read_parquet(RES / f"activity_{tag}_{what}_wide.parquet")


def activity_long() -> pd.DataFrame:
    return pd.read_parquet(RES / "activity_long.parquet")


def kinases() -> pd.DataFrame:
    return pd.read_csv(RES / "kinase_annotation.tsv", sep="\t")


def drug_targets() -> pd.DataFrame:
    return pd.read_csv(RES / "drug_targets.tsv", sep="\t")
