"""audit_F helpers: own condition classes, robust z, hand-picked canonical-marker positive controls (independent of PSP gold)."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load

CACHE = "/workspace/cache/audit_F"
LAB = {"26060313": "Cutillas", "31959955": "Cutillas", "32034161": "Olsen", "30459283": "Olsen", "36926954": "Kuster",
       "27794612": "Kuster", "28209619": "Kuster", "30901187": "Kuster", "29255136": "Ahn", "25850435": "Ahn",
       "25002506": "Roux", "27422710": "Haura", "34663829": "Mann", "38047086": "Trajanoski", "28188228": "Thibault",
       "29970458": "Heck", "26055452": "Burgess", "27134283": "Kettenbach", "30773488": "Sterneckert", "28165663": "Marto",
       "23898821": "Daub", "22115753": "Daub", "27732954": "Wu", "30459281": "Hallberg", "25404012": "Guha"}

CLASSES = {  # inhibitor class -> primary-target kinases (condition included if any primary target in set)
    "AKT": {"AKT1", "AKT2", "AKT3"}, "MEK": {"MAP2K1", "MAP2K2"}, "ERK": {"MAPK1", "MAPK3"}, "RAF": {"BRAF", "RAF1", "ARAF"},
    "RSK": {"RPS6KA1", "RPS6KA3"}, "CDK": {"CDK1", "CDK2", "CDK5", "CDK9"}, "PKC": {"PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCQ"},
    "MTOR": {"MTOR"}, "ALK": {"ALK", "ROS1"}, "EGFR": {"EGFR", "ERBB2"}, "PI3K": {"PIK3CA", "PIK3CB", "PIK3CD"},
}
# canonical readouts (not test sites); MAPK markers = ERK activation loop + ERK/RSK feedback/substrates
MARK = {
    "AKT": ["P49841_S9", "P49840_S21", "Q96B36_T246", "P49815_T1462", "O43524_S253", "Q12778_S256", "O60343_T642"],
    "MAPK": ["P27361_T202", "P27361_Y204", "P28482_T185", "P28482_Y187", "P04049_S296", "P04049_S289", "P04049_S301",
             "P15056_T401", "P15056_S750", "P15056_S151", "Q07889_S1134", "Q07889_S1161", "P00533_T693"],
    "RSKsub": ["Q07889_S1161", "P67809_S102", "P62753_S235", "P62753_S236"],
    "CDK": ["P06400_S807", "P06400_S811", "P06400_S780", "P06400_S795", "P06748_T199", "P02545_S22"],
    "PKC": ["P29966_S163", "P29966_S167", "P29966_S170", "Q15139_S742", "Q15139_S738"],
    "MTOR": ["Q13541_T37", "Q13541_T46", "Q13541_S65", "Q13541_T70", "O75385_S758", "Q96B36_S183", "P62753_S240", "P62753_S244"],
}
CLS2MARK = {"AKT": "AKT", "MEK": "MAPK", "ERK": "MAPK", "RAF": "MAPK", "RSK": "RSKsub", "CDK": "CDK", "PKC": "PKC", "MTOR": "MTOR",
            "ALK": "MAPK", "EGFR": "MAPK", "PI3K": "AKT"}


def conditions():
    C = load.conditions()
    C = C[C.qc_flag.eq("ok") & C.duplicate_of.isna() & C.perturbation_class.eq("kinase_inhibitor") & ~C.is_combination.astype(bool)].copy()
    pt = C.inhibitor_primary_targets.fillna("").str.split(";").apply(set)
    rows = []
    for cid, st, p, ag, rev, cl, t in zip(C.conditionID, C.study_pmid, pt, C.agents, C.contrast_reversed, C.cell_line_norm, C.time_min):
        if len(p) > 8:  # promiscuous
            continue
        for k, ks in CLASSES.items():
            if p & ks:
                rows.append((k, cid, st, LAB.get(st, "lab_" + st), ag, bool(rev), cl, t))
    return pd.DataFrame(rows, columns=["cls", "conditionID", "study", "lab", "agent", "rev", "cell", "time_min"])


def zmatrix(conds):
    p = f"{CACHE}/Z_inh.parquet"
    if os.path.exists(p):
        Z = pd.read_parquet(p)
        if set(conds) <= set(Z.columns):
            return Z
    X = load.sites_matrix(sorted(set(conds)))
    med = X.median(); mad = (X - med).abs().median() * 1.4826
    Z = ((X - med) / mad.clip(lower=0.05)).clip(-8, 8).astype(np.float32)
    Z.to_parquet(p)
    return Z


def marker_score(Z, CC):
    """Stouffer of class markers per condition (sign-corrected)."""
    out = {}
    for r in CC.itertuples():
        ms = MARK[CLS2MARK[r.cls]]
        v = Z[r.conditionID].reindex(ms).dropna() * (-1 if r.rev else 1)
        out[(r.cls, r.conditionID)] = (v.sum() / np.sqrt(len(v)) if len(v) >= 2 else np.nan, len(v))
    return out
