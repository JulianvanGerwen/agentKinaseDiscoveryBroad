"""Permutation-FDR calling of regulated dose-response curves.
Candidate curve: EC50 within tested range, |delta| >= log2(1.5), n>=6.  Statistic = F (4PL vs flat).
Per experiment FDR(t) = [#null candidates F>=t / NPERM] / #observed candidates F>=t ; call at q<=0.05 (monotone q).
Also a Hill test: F_hill = (sse_h1 - sse_inr)/(sse_inr/(n-4)) (Hill free vs Hill=1, both EC50 in range).
Output: results/b_I17/curves_all.parquet (all fits + q), results/b_I17/fdr_by_experiment.tsv
"""
import numpy as np, pandas as pd
OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"
F = pd.read_parquet(f"{OUT}/fits.parquet").reset_index(drop=True); N = pd.read_parquet(f"{OUT}/null_F.parquet")
E = pd.read_csv(f"{OUT}/experiments.tsv", sep="\t").set_index("exp_id")
THR = np.log2(1.5)
F["cand"] = F.ec_in_range & (F.delta.abs() >= THR) & (F.n >= 6) & (F.front.abs() < 1.5)
N["cand"] = N.ec_in_range & (N.delta.abs() >= THR) & (N.front.abs() < 1.5) if "front" in N else N.ec_in_range & (N.delta.abs() >= THR)
nperm = N.groupby("exp_id").perm.nunique()
F["q"] = 1.0
rows = []
for ex, g in F.groupby("exp_id"):
    nul = np.sort(N.loc[(N.exp_id == ex) & N.cand, "F"].values)
    c = g[g.cand]
    if len(c) == 0:
        continue
    t = c.F.values
    n_null = len(nul) - np.searchsorted(nul, t, side="left")
    rank_obs = np.array([(t >= x).sum() for x in t]) if len(t) < 2000 else (len(t) - np.searchsorted(np.sort(t), t, side="left"))
    fdr = (n_null / nperm[ex]) / np.maximum(rank_obs, 1)
    o = np.argsort(-t); q = np.minimum.accumulate(fdr[o][::-1])[::-1]; qq = np.empty_like(q); qq[o] = np.minimum(q, 1)
    F.loc[c.index, "q"] = qq
    rows.append(dict(exp_id=ex, n_sites=len(g), n_cand=len(c), n_null_cand_per_perm=len(nul) / nperm[ex],
                     n_sig=(qq <= 0.05).sum(), n_sig_down=((qq <= 0.05) & (c.delta < 0)).sum(), n_sig_up=((qq <= 0.05) & (c.delta > 0)).sum()))
F["sig"] = F.q <= 0.05
F["pEC50"] = 9 - F.logEC50
F["F_hill"] = (F.sse_h1 - F.sse_inr) / (F.sse_inr / np.maximum(F.n - 4, 1))
R = pd.DataFrame(rows).set_index("exp_id").join(E[["agent", "cell_line", "time_min", "pclass", "targets"]])
R["frac_sig"] = R.n_sig / R.n_sites
R.to_csv(f"{RES}/fdr_by_experiment.tsv", sep="\t")
F = F.join(E[["agent", "cell_line", "time_min", "pclass"]], on="exp_id")
F.to_parquet(f"{RES}/curves_all.parquet")
print(R.sort_values("frac_sig").to_string())
