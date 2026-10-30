"""Shared paths/helpers for I28 (phosphosite turnover atlas)."""
import numpy as np, pandas as pd
C28 = "/workspace/cache/b_I28"; R28 = "/workspace/results/b_I28"
C26 = "/workspace/cache/b_I26"
KSM = "<HACKATHON_ROOT>/kinaseSubMain"
FOCUS = ["AKT", "ERK", "MTOR", "CDK1/2/5/9", "EGFR/ERBB2", "ABL/KIT/PDGFR", "ALK/ROS1", "ROCK", "ATR", "ATM", "CK2", "PKA"]


def classify_units(L, fast_mean=-1.0, stable_mean=-0.5):
    """L: long (group, site, unit, study, z). Returns per (group,site) atlas stats."""
    L = L.assign(dec=L.z <= -1)
    st = L.groupby(["group", "site", "study"], observed=True).z.mean().groupby(["group", "site"], observed=True)
    A = L.groupby(["group", "site"], observed=True).agg(n_units=("z", "size"), n_studies=("study", "nunique"),
                                                        R=("z", "mean"), n_dec=("dec", "sum"), z_min=("z", "min")).reset_index()
    A = A.merge(st.mean().rename("R_study").reset_index(), on=["group", "site"])
    f = A.n_dec / A.n_units
    A["frac_dec"] = f
    A["cls"] = np.select([A.n_units < 2, (A.R <= fast_mean) & (A.n_dec >= 2) & (f >= 0.5), (A.R > stable_mean) & (f < 0.25)],
                         ["untestable", "fast", "stable"], "intermediate")
    return A
