"""Build dose-series experiment table (split merged series) and the Klaeger Kd matrix.
Outputs: cache/b_I17/experiments.tsv, cache/b_I17/klaeger_kd.parquet (protein x drug, nM; n.i. -> inf), results/b_I17/drug_map.tsv
"""
import os, re, sys
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load

OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"
C = load.conditions()
d = C[C.dose_series_id.notna()].copy()
d["idx"] = d.conditionID.str.extract(r"_(\d+)_log2$")[0].astype(float)
# experiment = dose_series_id + control_raw + (for decryptM) block split on repeated doses
rows = []
for ds, g in d.groupby("dose_series_id"):
    for ctrl, h in g.groupby("control_raw"):
        h = h.sort_values("idx")
        exp, seen, k = [], {}, 0
        # split when a dose 0 appears again after other doses (new sub-experiment)
        for _, r in h.iterrows():
            if r.dose_nM == 0 and len(exp) and any(x != 0 for x in seen.get(k, [])):
                k += 1
            seen.setdefault(k, []).append(r.dose_nM)
            exp.append(k)
        h = h.assign(block=exp)
        for b, hh in h.groupby("block"):
            nd = hh.dose_nM.nunique()
            rep = hh.dose_nM.value_counts().max()
            rows.append(dict(exp_id=f"{ds}_{re.sub('[^A-Za-z0-9]', '', str(ctrl))[:30]}_{b}", dose_series_id=ds,
                             study=hh.study_pmid.iloc[0], agent=hh.agents.iloc[0], cell_line=hh.cell_line_norm.iloc[0],
                             time_min=hh.time_min.iloc[0], pclass=hh.perturbation_class.iloc[0],
                             targets=hh.inhibitor_primary_targets.iloc[0], n_cond=len(hh), n_doses=nd, max_rep=rep,
                             conditions=";".join(hh.conditionID), doses=";".join(map(str, hh.dose_nM))))
E = pd.DataFrame(rows)
E.to_csv(f"{OUT}/experiments.tsv", sep="\t", index=False)
print(E.drop(columns=["conditions", "doses"]).to_string())

# Klaeger Kd matrix
x = pd.read_excel("/workspace/data/Klaeger/Table_S3.xlsx", sheet_name="Kinobeads Drugmatrix detailed")
x = x.set_index("Name")
meta = x[["Kinase", "Direct binder"]]
K = x.drop(columns=["Kinase", "Direct binder"])
def conv(v):
    if isinstance(v, str):
        v = v.strip()
        if v == "n.i.": return np.inf
        if v == "n.d.": return np.nan
        try: return float(v)
        except ValueError: return np.nan
    return float(v) if pd.notna(v) else np.nan
K = K.map(conv)
K.columns = [str(c) for c in K.columns]
K.to_parquet(f"{OUT}/klaeger_kd.parquet")
meta.to_csv(f"{OUT}/klaeger_meta.tsv", sep="\t")
norm = lambda s: re.sub("[^a-z0-9]", "", str(s).lower())
kn = {norm(c): c for c in K.columns}
alias = {"pd0325901": "pd325901", "azd8055": "azd8055", "mk2206": "mk2206"}
dm = []
for a in sorted(set(E.agent.dropna())):
    parts = a.split(";")
    m = [kn.get(alias.get(norm(p), norm(p))) for p in parts]
    dm.append(dict(agent=a, klaeger=";".join([str(v) for v in m]), all_mapped=all(v is not None for v in m)))
DM = pd.DataFrame(dm); DM.to_csv(f"{RES}/drug_map.tsv", sep="\t", index=False); print(DM.to_string())
print("Klaeger", K.shape, "finite", np.isfinite(K.values).sum())
