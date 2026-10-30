"""Phospho-response features per single-agent inhibitor condition.
Hijazi (31959955): centred within cell (site median across that cell's inhibitor conditions).
Other studies: uncentred (vs own control).  Out: results/b_I23/tables/T1_features.tsv"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/b_I23/tables/"; import os; os.makedirs(OUT, exist_ok=True)
C = load.conditions()
k = C[(C.perturbation_class == "kinase_inhibitor") & C.background_agents.isna() & (C.is_combination != True)
      & (C.qc_flag.fillna("ok") == "ok") & (C.contrast_reversed != True)].copy()
k = k[~k.kinase_inhibitors.str.contains(";", na=True)]
X = load.sites_matrix(k.conditionID)
EPS = 0.002
GROUPS = {"STE20_MAP4K": ["MAP4K1","MAP4K2","MAP4K3","MAP4K4","MAP4K5","TNIK","MINK1"],
          "NEK": [f"NEK{i}" for i in range(1, 12)], "IRAK4": ["IRAK4"],
          "MAP3K": None, "MARK_SIK": ["MARK1","MARK2","MARK3","MARK4","SIK1","SIK2","SIK3"],
          "AKT": ["AKT1","AKT2","AKT3"], "SFK": ["SRC","YES1","FYN","LYN","LCK","HCK","FGR","BLK"]}
Z = load.activity("corev3"); N = load.activity("corev3", "nsub")
# z_core wide
AL = load.activity_long(); AL = AL[(AL.method == "corev3") & AL.condition.isin(k.conditionID)]
ZC = AL.pivot(index="kinase", columns="condition", values="z_core"); NS = AL.pivot(index="kinase", columns="condition", values="n_substrates")
ZC = ZC.where(NS >= 3)
GROUPS["MAP3K"] = [g for g in ZC.index if g.startswith("MAP3K")]
D = load.drug_targets().set_index("agent")
motif = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")

def feats(x, thr=1.0):
    m = x.notna(); v = x[m]
    fu = (v >= thr).mean(); fd = (v <= -thr).mean()
    return dict(n_meas=int(m.sum()), f_up=fu, f_down=fd, LR=np.log((fu + EPS) / (fd + EPS)),
                f_up58=(v >= .58).mean(), f_down58=(v <= -.58).mean(),
                sum_up=v.clip(lower=0).sum() / len(v), sum_down=(-v).clip(lower=0).sum() / len(v))

rows = []
for (study, cell), g in k.groupby(["study_pmid", "cell_line_norm"]):
    Xs = X[g.conditionID]
    centred = study == "31959955"
    Xc = Xs.sub(Xs.median(axis=1), axis=0) if centred else Xs
    for _, r in g.iterrows():
        x = Xc[r.conditionID]
        f = feats(x)
        # random half-site control
        rng = np.random.default_rng(abs(hash(r.conditionID)) % 2**32)
        xs = x.dropna(); half = xs.sample(frac=.5, random_state=1)
        f["LR_half"] = np.log(((half >= 1).mean() + EPS) / ((half <= -1).mean() + EPS))
        # motif composition of up-layer: fraction of up sites with motif pct>=90 for group, minus same for all measured
        up = xs.index[xs >= 1]; allm = xs.index
        for gname, ks in GROUPS.items():
            ks2 = [q for q in ks if q in motif.columns]
            if ks2:
                mm = motif.reindex(allm)[ks2].max(axis=1) >= 90
                f[f"motif_{gname}"] = (mm.reindex(up).mean() - mm.mean()) if len(up) >= 10 else np.nan
            if r.conditionID in ZC.columns:
                zz = ZC.loc[ZC.index.intersection(ks), r.conditionID]
                f[f"z_{gname}"] = zz.mean() if zz.notna().any() else np.nan
        tg = str(D.primary_targets_in_universe.get(r.kinase_inhibitors, "")).split(";")
        if r.conditionID in ZC.columns:
            zt = ZC.loc[ZC.index.intersection(tg), r.conditionID]
            f["target_z"] = zt.mean() if zt.notna().any() else np.nan
        f.update(conditionID=r.conditionID, study=study, cell=cell, drug=r.kinase_inhibitors, dose_nM=r.dose_nM,
                 time_min=r.time_min, centred=centred, dose_series=r.dose_series_id)
        rows.append(f)
F = pd.DataFrame(rows)
F.to_csv(OUT + "T1_features.tsv", sep="\t", index=False)
print(F.shape); print(F[F.study == "31959955"][["f_up","f_down","LR","LR_half","target_z"]].describe())
