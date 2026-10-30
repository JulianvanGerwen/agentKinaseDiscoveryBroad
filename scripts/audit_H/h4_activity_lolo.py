"""AUDIT copy of I28 05_activity.py; MODE=loso (replicate) or lolo (turnover re-estimated without ALL studies of the evaluated study's lab).
I28 step 5: does turnover weighting of the PSP library improve activity inference on held-out inhibitor studies?
Scorer: PSP Stouffer of robust z (as infra psp_z / I26), >=3 substrates. Benchmark conditions = I26 set (inhibitor,
target confidence high/medium, qc ok, dose>=100nM). Turnover is re-estimated WITHOUT the evaluated study (LOSO).
Variants (edges of kinases in atlas families only; other edges untouched):
  base; drop_stable (cls==stable removed); keep_fast_only (measured non-fast removed);
  weighted (w = clip(-R, 0.25, 4) for measured edges, 1 otherwise; weighted Stouffer sum(wz)/sqrt(sum w^2));
  rand_drop (drop same number of measured edges per kinase at random, 10 reps); rand_weight (weights permuted within kinase, 10 reps)
Out: tables/activity_turnover_conditions.tsv, tables/activity_turnover_summary.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import rankdata, wilcoxon
sys.path.insert(0, "/workspace/scripts/b_I28"); from lib28 import *
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
rng = np.random.default_rng(5)
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
lib = lib[lib.site.isin(Z.index)][["kinase", "site"]].drop_duplicates()
G = pd.read_csv(f"{C26}/groups.tsv", sep="\t"); k2g = dict(zip(G.kinase, G.group))
lib["group"] = lib.kinase.map(k2g)
L = pd.read_parquet(f"{C28}/L_units.parquet"); L = L[L.set == "substrate"]
L = L[L.group != "MEK"]  # ERK pathway handled via ERK; MEK/RAF libraries are tiny
sites = sorted(lib.site.unique()); sidx = {s: i for i, s in enumerate(sites)}
Zs = Z.reindex(sites); fin = Zs.notna().values; Zv = np.where(fin, Zs.values, 0.0)
kins = sorted(lib.kinase.unique()); kidx = {k: i for i, k in enumerate(kins)}
conds = list(Zs.columns); cinfo = sel.set_index("conditionID")
tg = {c: set(cinfo.at[c, "inhibitor_primary_targets_in_universe"].split(";")) for c in conds}
KI = np.array([kidx[k] for k in lib.kinase]); SI = np.array([sidx[s] for s in lib.site])


def score(w):
    K = np.zeros((len(kins), len(sites)), np.float32); K[KI, SI] = w
    N = (K > 0).astype(np.float32) @ fin; S = K @ Zv; Q = (K ** 2) @ fin
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(N >= 3, S / np.sqrt(Q), np.nan)


def auroc(z, t):
    ok = ~np.isnan(z); z, t = z[ok], t[ok]
    nt, nn = t.sum(), (~t).sum()
    if nt == 0 or nn == 0: return np.nan
    r = rankdata(-z)
    return (r[t].sum() - nt * (nt + 1) / 2) / (nt * nn)


def variants(A):
    m = lib.merge(A[["group", "site", "R", "cls"]], on=["group", "site"], how="left")
    meas = m.cls.isin(["fast", "intermediate", "stable"]).values
    out = {"base": np.ones(len(m))}
    out["drop_stable"] = np.where(m.cls.eq("stable"), 0, 1.0)
    out["keep_fast_only"] = np.where(meas & ~m.cls.eq("fast"), 0, 1.0)
    wt = np.where(meas, np.clip(-m.R.fillna(-1).values, 0.25, 4), 1.0); out["weighted"] = wt
    for r in range(10):
        d = np.ones(len(m)); dk = np.ones(len(m)); w = wt.copy()
        for k, idx in m[meas].groupby("kinase").groups.items():
            idx = np.array(list(idx)); nst = int(m.loc[idx, "cls"].eq("stable").sum())
            if nst: d[rng.choice(idx, nst, replace=False)] = 0
            nnf = int((~m.loc[idx, "cls"].eq("fast")).sum())
            if nnf: dk[rng.choice(idx, nnf, replace=False)] = 0
            w[idx] = rng.permutation(wt[idx])
        out[f"rand_drop_{r}"] = d; out[f"rand_keepn_{r}"] = dk; out[f"rand_weight_{r}"] = w
    edited = set(m.kinase[meas])
    return out, edited


MODE = os.environ.get("MODE", "lolo")
LAB = {"31959955": "Cutillas", "26060313": "Cutillas", "36926954": "Kuster", "27794612": "Kuster", "28209619": "Kuster",
       "30901187": "Kuster", "25850435": "Ahn", "29255136": "Ahn", "30459283": "Olsen", "32034161": "Olsen",
       "22115753": "Daub", "23898821": "Daub", "21712546": "Kettenbach", "27134283": "Kettenbach", "29214152": "Kettenbach",
       "25003641": "Tichy", "30001349": "Tichy"}
lab_of = lambda x: LAB.get(x, "L" + x)
L["lab"] = L.study.map(lab_of)
R28 = "/workspace/results/audit_H"
rows = []
studies = sorted(sel.study_pmid.unique())
Afull = classify_units(L)
for st in studies:
    excl = (L.study != st) if MODE == "loso" else (L.lab != lab_of(st))
    A = classify_units(L[excl])
    A = A[A.n_units >= 2]
    vs, edited = variants(A)
    cs = [c for c in sel.conditionID[sel.study_pmid == st]]
    ci = [conds.index(c) for c in cs]
    S = {v: score(w)[:, ci] for v, w in vs.items()}
    for j, c in enumerate(cs):
        t = np.array([k in tg[c] for k in kins])
        both = np.all([~np.isnan(S[v][:, j]) for v in S], axis=0)
        has_edit = bool(tg[c] & edited)
        for v in S:
            rows.append((st, c, v, auroc(np.where(both, S[v][:, j], np.nan), t), has_edit))
    print(st, len(cs), flush=True)
R = pd.DataFrame(rows, columns=["study", "condition", "variant", "auroc", "target_edited"])
R["vbase"] = R.variant.str.replace(r"_\d+$", "", regex=True)
R.to_csv(f"{R28}/h4_activity_conditions_{MODE}.tsv", sep="\t", index=False)
P = R.groupby(["study", "condition", "vbase", "target_edited"]).auroc.mean().unstack("vbase").reset_index()
out = []
for sub, d in [("all", P), ("target_edited", P[P.target_edited])]:
    d = d.dropna(subset=["base"])
    for v in ["drop_stable", "keep_fast_only", "weighted", "rand_drop", "rand_weight", "rand_keepn"]:
        delta = d[v] - d.base
        ds = (d.assign(delta=delta).groupby("study").delta.mean()).dropna()
        bs = [rng.choice(ds.values, len(ds)).mean() for _ in range(2000)]
        out.append(dict(subset=sub, variant=v, n_cond=delta.notna().sum(), n_studies=len(ds), base_auroc=d.base.mean(),
                        mean_delta_cond=delta.mean(), mean_delta_study=ds.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                        wilcoxon_p_study=wilcoxon(ds).pvalue if (ds != 0).sum() > 5 else np.nan))
O = pd.DataFrame(out); O.to_csv(f"{R28}/h4_activity_summary_{MODE}.tsv", sep="\t", index=False)
print(O.round(4).to_string())
