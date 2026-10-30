"""I03 step 8: per-edge evidence detail for curated candidates: per-study drop (agents, n conditions, mean z, mean log2FC),
tumour r per cohort, dose-series hits, and family/kinase-swap comparison at the same site.
Input: list of (kinase, site) in CANDS. Output results/b_I03/tables/top_edges_evidence_detail.tsv and top_edges_kinase_swap.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/cache/b_I03"; RES = "/workspace/results/b_I03/tables"
CANDS = [l.split() for l in open(f"{OUT}/curated_candidates.txt") if l.strip() and not l.startswith("#")]
KC = pd.read_csv(f"{OUT}/kinase_conditions.tsv", sep="\t")
C = load.conditions().set_index("conditionID")
sites = sorted({s for _, s in CANDS})
X = load.sites_matrix(sorted(set(KC.conditionID)), sites=sites)
med = load.sites_matrix(sorted(set(KC.conditionID))).median() if False else None
# robust z using the full columns (recompute med/MAD per condition)
full = load.sites_matrix(sorted(set(KC.conditionID)))
md = full.median(); mad = (full - md).abs().median() * 1.4826
Z = (X - md) / mad.clip(lower=0.05)
del full
rows = []
for k, s in CANDS:
    kc = KC[KC.kinase == k]
    for st, d in kc.groupby("study"):
        cc = [c for c in d.conditionID if c in Z.columns and pd.notna(Z.loc[s, c])] if s in Z.index else []
        if not cc: continue
        sg = d.set_index("conditionID").sign.reindex(cc).values
        z = Z.loc[s, cc].values * sg; fc = X.loc[s, cc].values * sg
        rows.append(dict(kinase=k, site=s, study=st, agents=";".join(sorted(set(d.set_index("conditionID").agent.reindex(cc)))),
                         cells=";".join(sorted(set(C.cell_line_norm.reindex(cc).dropna().astype(str)))), n_cond=len(cc),
                         mean_z=np.mean(z), mean_log2fc=np.mean(fc), min_log2fc=np.min(fc)))
D = pd.DataFrame(rows); D.to_csv(f"{RES}/top_edges_evidence_detail.tsv", sep="\t", index=False, float_format="%.3g")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(D.round(2).to_string())
# kinase swap at the same site: all panel kinases with >=2 studies
P = pd.read_parquet(f"{OUT}/scored.parquet")
P = P[P.site.isin(sites)][["site", "kinase", "family", "n_st", "n_st_drop", "stouffer", "ks_dev_pct", "motif", "tum_z", "dose_r"]]
P.sort_values(["site", "stouffer"]).to_csv(f"{RES}/top_edges_kinase_swap.tsv", sep="\t", index=False, float_format="%.3g")
