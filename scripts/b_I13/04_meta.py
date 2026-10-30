"""Meta-analysis of stage-2 coefficients with empirical (random-gene CNV) nulls."""
import numpy as np, pandas as pd
from scipy import stats
import sys
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
R = "/workspace/results/b_I13"
S2 = pd.read_parquet(f"{R}/cache/stage2_all{TAG}.parquet")
S2["kind"] = S2.exposure.str.split(":").str[0].where(S2.exposure.str.startswith("null_"), "real")
S2["nulltype"] = np.where(S2.exposure.str.contains("_prot"), "null_prot", np.where(S2.exposure.str.contains("hotspot"), "null_perm", "null_other"))
def dl(b, se):
    w = 1 / se ** 2; bf = (w * b).sum() / w.sum(); q = (w * (b - bf) ** 2).sum(); k = len(b)
    tau2 = max(0, (q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    ws = 1 / (se ** 2 + tau2); br = (ws * b).sum() / ws.sum(); ser = np.sqrt(1 / ws.sum())
    return br, ser, tau2, q
rows, per = [], []
nul = S2[S2.kind.str.startswith("null") & (S2.kind != "null_8p")].groupby(["cohort", "test", "kind"]).coef.agg(["mean", "std"]).rename(columns={"mean": "null_mean", "std": "null_sd"}).reset_index().rename(columns={"kind": "nulltype"})
nul8 = S2[S2.kind == "null_8p"].groupby(["cohort", "test"]).coef.agg(["mean", "std"]).rename(columns={"mean": "null8p_mean", "std": "null8p_sd"})
R2 = S2[S2.kind == "real"].merge(nul, on=["cohort", "test", "nulltype"], how="left").join(nul8, on=["cohort", "test"])
R2["z_emp"] = (R2.coef - R2.null_mean) / R2.null_sd
R2["z_vs8p"] = (R2.coef - R2.null8p_mean) / R2.null8p_sd
R2.to_csv(f"{R}/tables/t2_stage2_per_cohort{TAG}.tsv", sep="\t", index=False)
# meta: DL random effects on (coef - null_mean) with SE = null_sd (empirical); also Stouffer on z_emp
# null calibration of the meta: apply same meta to each null gene index position across cohorts (random pairing)
for (test, ex), d in R2.groupby(["test", "exposure"]):
    d = d.dropna(subset=["null_sd"])
    if len(d) == 0: continue
    b, se, tau2, q = dl((d.coef - d.null_mean).values, d.null_sd.values)
    zs = d.z_emp.sum() / np.sqrt(len(d))
    rows.append(dict(test=test, exposure=ex, k=len(d), beta_re=b, ci_lo=b - 1.96 * se, ci_hi=b + 1.96 * se, z_re=b / se, p_re=2 * stats.norm.sf(abs(b / se)),
                     tau2=tau2, stouffer_z=zs, n_cohorts_pos=int((d.z_emp > 0).sum()), n_cohorts_sig=int((d.z_emp.abs() > 1.96).sum()),
                     mean_z_vs8p=d.z_vs8p.mean(), cohorts=";".join(f"{c}:{z:.1f}" for c, z in zip(d.cohort, d.z_emp))))
M = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
M["fdr_bh"] = multipletests(M.p_re, method="fdr_bh")[1]
M.sort_values("p_re").to_csv(f"{R}/tables/t2_meta{TAG}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
PRE = {"PPP2R2A_cnv": ["B55_lit|ST", "bipartite_basic|ST", "isT|ST", "P1|ST"], "PPP2R2A_prot": ["B55_lit|ST", "bipartite_basic|ST", "isT|ST", "P1|ST"], "PPP2": ["bipartite_basic|ST", "B55_lit|ST", "B56_lit|ST"], "ARPP|ENSA|MASTL|PPP1": ["bipartite_basic|ST", "B55_lit|ST"], "XPPP2R2A": ["B55_lit|ST", "isT|ST", "P1|ST", "bipartite_basic|ST", "basic_m3p3|ST", "acid_m3p3|ST"], "PPP2R1A": ["B56_lit|ST", "R_m3|ST", "B55_lit|ST", "isT|ST"],
       "PPM1D": ["Q1|ST", "SQ_acid|ST"], "PTP": ["isY|all"], "PTEN": ["akt_motif|ST"], "PPP6C": ["isT|ST", "B55_lit|ST"]}
for g, tests in PRE.items():
    sel = M[M.exposure.str.match(g) & M.test.isin(tests)] if g != "PTP" else M[M.exposure.str.match("PTP(RD|N1|N11|N12|N2|N6)") & M.test.isin(tests)]
    print(sel[["exposure", "test", "k", "beta_re", "ci_lo", "ci_hi", "p_re", "fdr_bh", "stouffer_z", "n_cohorts_pos", "mean_z_vs8p", "cohorts"]].round(4).to_string(index=False))
