"""Rank real protein exposures against 1000 random proteins measured in all 10 cohorts (same protein across cohorts = strict null)."""
import numpy as np, pandas as pd
import sys
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
R = "/workspace/results/b_I13"
S2 = pd.read_parquet(f"{R}/cache/stage2_all{TAG}.parquet")
np_ = S2[S2.exposure.str.startswith("null_prot")]
st = np_.groupby(["cohort", "test"]).coef.agg(["mean", "std"])
d = S2[S2.exposure.str.startswith("null_prot") | (S2.exposure.str.contains("_prot") & ~S2.exposure.str.startswith("null"))].join(st, on=["cohort", "test"])
d["z"] = (d.coef - d["mean"]) / d["std"]
M = d.groupby(["test", "exposure"]).agg(k=("z", "size"), stouffer=("z", lambda z: z.sum() / np.sqrt(len(z))), npos=("z", lambda z: int((z > 0).sum())), mean_coef=("coef", "mean")).reset_index()
M = M[M.k >= 8]
out = []
for test, g in M.groupby("test"):
    nul = g[g.exposure.str.startswith("null_prot")].stouffer.values
    for _, r in g[~g.exposure.str.startswith("null")].iterrows():
        out.append(dict(test=test, exposure=r.exposure, k=r.k, stouffer=r.stouffer, npos=r.npos, n_null=len(nul),
                        emp_p_two_sided=(1 + (np.abs(nul) >= abs(r.stouffer)).sum()) / (1 + len(nul)),
                        rank_pct=(nul < r.stouffer).mean(), null_sd=nul.std()))
O = pd.DataFrame(out).sort_values("emp_p_two_sided"); O.to_csv(f"{R}/tables/t2_protein_exposure_rank_vs_1000_proteins{TAG}.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(O[O.test.isin(["bipartite_basic|ST", "B55_lit|ST", "isT|ST", "P1|ST", "isY|all", "B56_lit|ST", "basic_m3p3|ST"])].head(30).round(4).to_string(index=False))
print(O[O.exposure.str.startswith("PPP2R2A")].round(4).to_string(index=False))
# also which null proteins score top for bipartite_basic
g = M[(M.test == "bipartite_basic|ST") & M.exposure.str.startswith("null_prot")].sort_values("stouffer", ascending=False)
print(g.head(15).to_string(index=False))
