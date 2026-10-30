"""I26 step 5: consequence for activity inference (PSP mean-z / Stouffer, same as infra psp_z).
Benchmark: infra inhibitor-target sanity conditions (kinase_inhibitor, subclass inhibitor / +stimulus / reversed, target
confidence high/medium, qc ok, dose>=100 nM or unknown). Per condition AUROC of curated primary targets among scored kinases.
Edge audit is re-derived WITHOUT the evaluated condition's study (leave-one-study-out, LOSO) and, separately, in 20 random
half splits (audit half -> evaluate other half). Variants:
  base      PSP library (infra psp_ks_edges_human_mapped)
  prune     drop edges whose (group,site) is 'unsupported'
  prune+amb drop 'unsupported' and 'ambiguous'
  prune+add prune + add verified non-upstream alternative edges (to all members of the alternative group)
  rand      drop the same number of edges per kinase at random among that kinase's measured, non-supported audited edges (20 reps)
  rand_any  drop the same number at random among all of that kinase's library edges (20 reps)
min substrates = 3 for all variants; AUROC computed on kinases scored in base and the variant.
Outputs tables/activity_benchmark_conditions.tsv, tables/activity_benchmark_summary.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import rankdata, wilcoxon
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load

rng = np.random.default_rng(7)
C = load.conditions()
sel = C[C.perturbation_class.eq("kinase_inhibitor")
        & C.perturbation_subclass.isin(["inhibitor", "inhibitor+stimulus", "inhibitor_reversed_contrast"])
        & C.target_confidence.isin(["high", "medium"]) & C.qc_flag.eq("ok")
        & (C.dose_nM.isna() | C.dose_nM.ge(100)) & C.inhibitor_primary_targets_in_universe.notna()].copy()
sign = np.where(sel.perturbation_subclass.eq("inhibitor_reversed_contrast"), -1.0, 1.0)
X = load.sites_matrix(sel.conditionID)
med = X.median(axis=0); mad = (X - med).abs().median(axis=0) * 1.4826
Z = ((X - med) / mad) * sign
lib = pd.read_csv("/workspace/results/infra_compendium/psp_ks_edges_human_mapped.tsv", sep="\t")
lib = lib[lib.site.isin(Z.index)]
G = pd.read_csv(f"{C26}/groups.tsv", sep="\t"); k2g = dict(zip(G.kinase, G.group)); members = G.groupby("group").kinase.apply(list).to_dict()
lib["group"] = lib.kinase.map(k2g)
Lm = pd.read_parquet(f"{C26}/L_main.parquet"); Lm["group"] = Lm.group.astype(str); Lm["study"] = Lm.study.astype(str)
Lm_psp = Lm.merge(lib[["group", "site"]].drop_duplicates(), on=["group", "site"])
CORR = pd.read_csv(f"{R26}/tables/corrections.tsv", sep="\t")
UP = {"EGFR/ERBB2", "MET", "ALK/ROS1", "ABL/KIT/PDGFR", "RAF", "MEK"}
sites = sorted(lib.site.unique()); sidx = {s: i for i, s in enumerate(sites)}
Zs = Z.reindex(sites)
fin = Zs.notna().values; Zv = np.where(fin, Zs.values, 0.0)
kins = sorted(set(lib.kinase) | set(G.kinase)); kidx = {k: i for i, k in enumerate(kins)}


def score(edges):
    K = np.zeros((len(kins), len(sites)), np.float32)
    e = edges[edges.site.isin(sidx)]
    K[[kidx[k] for k in e.kinase], [sidx[s] for s in e.site]] = 1
    N = K @ fin; S = K @ Zv
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(N >= 3, S / np.sqrt(N), np.nan)       # kinases x conditions


def auroc(z, is_t):
    ok = ~np.isnan(z)
    z, is_t = z[ok], is_t[ok]
    nt, nn = is_t.sum(), (~is_t).sum()
    if nt == 0 or nn == 0:
        return np.nan, np.nan
    r = rankdata(-z)
    return (r[is_t].sum() - nt * (nt + 1) / 2) / (nt * nn), z[is_t].mean()


def variants_for(A):
    """A = classification of PSP (group,site) from audit studies -> dict of edge libraries."""
    cl = lib.merge(A[["group", "site", "cls"]], on=["group", "site"], how="left")
    out = {"base": lib[["kinase", "site"]]}
    out["prune"] = cl[cl.cls != "unsupported"][["kinase", "site"]]
    out["prune+amb"] = cl[~cl.cls.isin(["unsupported", "ambiguous"])][["kinase", "site"]]
    un = cl[cl.cls == "unsupported"][["group", "site"]].drop_duplicates()
    # corrections: alt group supported (verified) for unsupported site, recomputed from same audit A
    sup = A[A.cls == "supported"][["group", "site"]].rename(columns={"group": "alt"})
    add = un.merge(sup, on="site"); add = add[(add.alt != add.group) & ~add.alt.isin(UP)]
    addk = pd.DataFrame([(k, s) for a, s in zip(add.alt, add.site) for k in members[a]], columns=["kinase", "site"])
    out["prune+add"] = pd.concat([out["prune"], addk]).drop_duplicates()
    nrem = cl[cl.cls == "unsupported"].groupby("kinase").size()
    pool_t = cl[cl.cls.isin(["unsupported", "ambiguous"])]
    for rep in range(10):
        drop = []
        for k, n in nrem.items():
            p = pool_t[pool_t.kinase == k].site.values
            p_any = lib[lib.kinase == k].site.values
            drop += [("rand", k, s) for s in rng.choice(p, n, replace=False)]
            drop += [("rand_any", k, s) for s in rng.choice(p_any, n, replace=False)]
        D = pd.DataFrame(drop, columns=["v", "kinase", "site"])
        for v in ["rand", "rand_any"]:
            d = D[D.v == v]
            out[f"{v}_{rep}"] = lib[~(lib.kinase + "|" + lib.site).isin(set(d.kinase + "|" + d.site))][["kinase", "site"]]
    return out, len(un)


conds = list(Zs.columns)
cinfo = sel.set_index("conditionID")
tg = {c: set(cinfo.at[c, "inhibitor_primary_targets_in_universe"].split(";")) for c in conds}


def evaluate(edgesets, cond_subset, tag):
    res = []
    ci = [conds.index(c) for c in cond_subset]
    scores = {v: score(e)[:, ci] for v, e in edgesets.items()}
    for j, c in enumerate(cond_subset):
        t = np.array([k in tg[c] for k in kins])
        both = ~np.isnan(scores["base"][:, j])
        for v in edgesets:
            both &= ~np.isnan(scores[v][:, j])
        for v in edgesets:
            z = np.where(both, scores[v][:, j], np.nan)
            a, mz = auroc(z, t)
            res.append((tag, c, v, a, mz))
    return res


rows = []
# ---- LOSO
for st, d in sel.groupby("study_pmid"):
    A = classify(Lm_psp, studies_excl=[st])
    es, nun = variants_for(A)
    rows += evaluate(es, list(d.conditionID), "LOSO")
    print("LOSO", st, len(d), "unsupported units", nun, flush=True)
# ---- half splits
studies = sorted(Lm.study.unique())
for sp in range(20):
    perm = rng.permutation(studies); half = set(perm[: len(perm) // 2])
    A = classify(Lm_psp, studies_only=half)
    es, nun = variants_for(A)
    ev = sel[~sel.study_pmid.isin(half)]
    rows += evaluate(es, list(ev.conditionID), f"split{sp}")
    print("split", sp, "audit studies", len(half), "unsupported units", nun, "eval conds", len(ev), flush=True)
R = pd.DataFrame(rows, columns=["scheme", "condition", "variant", "auroc", "mean_target_z"])
R["study"] = R.condition.map(cinfo.study_pmid); R["target_confidence"] = R.condition.map(cinfo.target_confidence)
R["targets"] = R.condition.map(cinfo.inhibitor_primary_targets_in_universe)
R.to_parquet(f"{C26}/activity_benchmark_long.parquet")
