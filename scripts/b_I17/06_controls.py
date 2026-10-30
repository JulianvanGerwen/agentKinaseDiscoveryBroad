"""Controls and replication.
1. Replicate / cross-cell-line agreement of pEC50 and Hill for same drug (sites significant in both).
2. Drug-label shuffle: experiment-level Spearman(pKd_target, median pEC50 of down sites); site-level median |shift|.
3. Hill-slope calibration: simulate every significant curve under Hill=1 (same EC50, amplitude, residual SD,
   dose design) and refit -> distribution of estimated Hill under the null of no ultrasensitivity.
Outputs: results/b_I17/replication_pairs.tsv, controls_shuffle.tsv, hill_calibration.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, "/workspace/scripts/b_I17"); from fitlib import fit_matrix, HILLS
OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"
F = pd.read_parquet(f"{RES}/curves_all.parquet")
E = pd.read_csv(f"{OUT}/experiments.tsv", sep="\t").set_index("exp_id")
A = pd.read_csv(f"{RES}/shift_sites.tsv", sep="\t")
X = pd.read_csv(f"{RES}/shift_by_experiment.tsv", sep="\t")
rng = np.random.default_rng(3)

# 1. replication
ki = [ex for ex in set(F.exp_id) if E.loc[ex, "pclass"] == "kinase_inhibitor"]
pairs = []
for i, a in enumerate(sorted(ki)):
    for b in sorted(ki)[i + 1:]:
        if E.loc[a, "agent"] != E.loc[b, "agent"]: continue
        fa = F[(F.exp_id == a)].set_index("site"); fb = F[(F.exp_id == b)].set_index("site")
        both = fa.index.intersection(fb.index)
        sa, sb = fa.loc[both], fb.loc[both]
        s2 = both[sa.sig.values & sb.sig.values & (np.sign(sa.delta.values) == np.sign(sb.delta.values))]
        s_any = both[sa.sig.values | sb.sig.values]
        kind = "replicate" if (E.loc[a, "cell_line"] == E.loc[b, "cell_line"] and E.loc[a, "time_min"] == E.loc[b, "time_min"]) else \
               "time" if E.loc[a, "cell_line"] == E.loc[b, "cell_line"] else "cell_line"
        r = dict(agent=E.loc[a, "agent"], exp_a=a, exp_b=b, kind=kind, cl_a=E.loc[a, "cell_line"], cl_b=E.loc[b, "cell_line"],
                 n_common=len(both), n_sig_a=int(sa.sig.sum()), n_sig_b=int(sb.sig.sum()), n_sig_both=len(s2),
                 jaccard=len(s2) / max(len(s_any), 1))
        if len(s2) >= 10:
            pa, pb = fa.loc[s2, "pEC50"], fb.loc[s2, "pEC50"]
            ha, hb = np.log2(fa.loc[s2, "hill"]), np.log2(fb.loc[s2, "hill"])
            r.update(rho_pEC50=spearmanr(pa, pb)[0], med_abs_dpEC50=(pa - pb).abs().median(), frac_within_0p5=((pa - pb).abs() <= 0.5).mean(),
                     rho_log2hill=spearmanr(ha, hb)[0],
                     steep_agree=((fa.loc[s2, "hill"] >= 2) == (fb.loc[s2, "hill"] >= 2)).mean(),
                     frac_steep_a=(fa.loc[s2, "hill"] >= 2).mean())
        pairs.append(r)
Pr = pd.DataFrame(pairs); Pr.to_csv(f"{RES}/replication_pairs.tsv", sep="\t", index=False)
print(Pr.round(3).to_string())

# 2. drug-label shuffle
rows = []
Xd = X.dropna(subset=["median_pEC50_down"])
drug_pkd = Xd.groupby("agent").pKd_target.first()
obs = spearmanr(Xd.pKd_target, Xd.median_pEC50_down)[0]
dn = A[A.direction == "down"]
obs_abs = dn["shift"].abs().median()
null_r, null_abs = [], []
agents = drug_pkd.index.values
for p in range(2000):
    perm = dict(zip(agents, rng.permutation(drug_pkd.values)))
    pk = Xd.agent.map(perm)
    null_r.append(spearmanr(pk, Xd.median_pEC50_down)[0])
    if p < 500:
        null_abs.append((dn.pEC50 - dn.agent.map(perm)).abs().median())
null_r = np.array(null_r); null_abs = np.array(null_abs)
rows.append(dict(test="exp-level Spearman(pKd_target, median pEC50 down)", observed=obs, null_mean=null_r.mean(),
                 p_perm=(np.sum(null_r >= obs) + 1) / (len(null_r) + 1), n=len(Xd), n_drugs=Xd.agent.nunique()))
rows.append(dict(test="site-level median |shift| (down)", observed=obs_abs, null_mean=null_abs.mean(),
                 p_perm=(np.sum(null_abs <= obs_abs) + 1) / (len(null_abs) + 1), n=len(dn), n_drugs=dn.agent.nunique()))
# drug-level (one point per drug) version
dd = Xd.groupby("agent").agg(pKd=("pKd_target", "first"), med=("median_pEC50_down", "median"))
rd = spearmanr(dd.pKd, dd.med)[0]
nd = np.array([spearmanr(rng.permutation(dd.pKd.values), dd.med)[0] for _ in range(5000)])
rows.append(dict(test="drug-level Spearman(pKd_target, median pEC50 down)", observed=rd, null_mean=nd.mean(),
                 p_perm=(np.sum(nd >= rd) + 1) / (len(nd) + 1), n=len(dd), n_drugs=len(dd)))
C = pd.DataFrame(rows); C.to_csv(f"{RES}/controls_shuffle.tsv", sep="\t", index=False); print(C.to_string())

# 3. Hill calibration under Hill=1
cal = []
for ex, g in A.groupby("exp_id"):
    doses = np.array([float(v) for v in E.loc[ex, "doses"].split(";")])
    nz = doses[doses > 0]; lo, hi = np.log10(nz.min()), np.log10(nz.max())
    lx = np.where(doses > 0, np.log10(np.maximum(doses, 1e-12)), lo - 1.0)
    sd = np.sqrt(g.sse.values / np.maximum(g.n.values - 4, 1))
    for rep in range(3):
        s = 1 / (1 + 10 ** (1.0 * (g.logEC50.values[:, None] - lx[None, :])))
        Y = g.front.values[:, None] + g.delta.values[:, None] * s + rng.normal(size=(len(g), len(lx))) * sd[:, None]
        r = fit_matrix(Y, lx, lo, hi)
        cal.append(pd.DataFrame(dict(exp_id=ex, site=g.site.values, hill_obs=g.hill.values, hill_sim_h1=r["hill"],
                                     F_hill_sim=(r["sse_h1"] - r["sse_inr"]) / (r["sse_inr"] / np.maximum(g.n.values - 4, 1)))))
Cal = pd.concat(cal); Cal.to_parquet(f"{RES}/hill_calibration.parquet")
s = dict(frac_obs_hill_ge2=(A.hill >= 2).mean(), frac_sim_h1_hill_ge2=(Cal.hill_sim_h1 >= 2).mean(),
         frac_obs_hill_le0p7=(A.hill <= 0.7).mean(), frac_sim_h1_hill_le0p7=(Cal.hill_sim_h1 <= 0.7).mean(),
         median_obs_hill=A.hill.median(), median_sim_hill=Cal.hill_sim_h1.median(),
         Fhill_q95_sim=np.quantile(Cal.F_hill_sim, .95), frac_obs_Fhill_gt_simq95=(A.F_hill > np.quantile(Cal.F_hill_sim, .95)).mean())
pd.Series(s).to_csv(f"{RES}/hill_calibration.tsv", sep="\t"); print(pd.Series(s))
