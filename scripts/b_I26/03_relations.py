"""I26 step 3: per-kinase support rates and relation to PSP evidence type, references, family and motif.
Unit = (inhibitor group, site) since evidence is family-level. Outputs tables/per_group_support.tsv,
tables/per_kinase_support.tsv, tables/relations_models.tsv, tables/relations_strata.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, mannwhitneyu
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *

EA = pd.read_csv(f"{R26}/tables/edge_audit.tsv", sep="\t")
T = EA[EA.cls != "untestable"]
U = T.groupby(["group", "site"]).agg(cls=("cls", "first"), mean_z=("mean_z", "first"), n_meas=("n_meas", "first"),
                                     n_dec=("n_dec", "first"), in_vivo=("in_vivo", "max"), in_vitro=("in_vitro", "max"),
                                     LT=("LT_LIT", "max"), MS=("MS_LIT", "max"), n_records=("n_records", "sum"),
                                     n_reg_pmids=("n_reg_pmids", "max"), motif=("motif_pct", "max"),
                                     mgroup=("manning_group", "first"), kinases=("kinase", lambda x: ";".join(sorted(x))),
                                     gene=("gene", "first"), is_gold=("is_gold", "max")).reset_index()
U["supported"] = (U.cls == "supported").astype(int)
U["evidence"] = np.select([U.in_vivo & U.in_vitro, U.in_vivo], ["vivo+vitro", "vivo_only"], "vitro_only")
U["logLT"] = np.log1p(U.LT); U["logMS"] = np.log1p(U.MS)
U["motif_bin"] = pd.cut(U.motif, [-1, 50, 80, 90, 101], labels=["<50", "50-80", "80-90", ">=90"])
U["LT_bin"] = pd.cut(U.LT, [-1, 0, 2, 9, 1e6], labels=["0 (HTP only)", "1-2", "3-9", ">=10"])
U.to_csv(f"{C26}/units_testable.tsv", sep="\t", index=False)

# per group / per kinase rates
rc = pd.read_csv(f"{R26}/tables/random_control.tsv", sep="\t").set_index("group")
pg = U.groupby("group").agg(n_testable=("cls", "size"), n_supported=("supported", "sum"),
                            n_unsupported=("cls", lambda x: (x == "unsupported").sum()),
                            n_ambiguous=("cls", lambda x: (x == "ambiguous").sum()), median_mean_z=("mean_z", "median"))
pg["support_rate"] = pg.n_supported / pg.n_testable; pg["unsupported_rate"] = pg.n_unsupported / pg.n_testable
pg["rand_support_rate"] = rc.n_rand_exp_supported / rc.n_edges_matched
pg["rand_unsupported_rate"] = rc.n_rand_exp_unsupported / rc.n_edges_matched
pg["n_psp_edges_total"] = EA.groupby("group").size()
pg.sort_values("n_testable", ascending=False).to_csv(f"{R26}/tables/per_group_support.tsv", sep="\t", float_format="%.3f")
print(pg.sort_values("n_testable", ascending=False).round(3).to_string())
pk = T.groupby("kinase").agg(group=("group", "first"), n_testable=("cls", "size"), supported=("cls", lambda x: (x == "supported").sum()),
                             unsupported=("cls", lambda x: (x == "unsupported").sum()), ambiguous=("cls", lambda x: (x == "ambiguous").sum()))
pk["support_rate"] = pk.supported / pk.n_testable; pk["unsupported_rate"] = pk.unsupported / pk.n_testable
pk["n_psp_edges_total"] = EA.groupby("kinase").size()
pk.sort_values("n_testable", ascending=False).to_csv(f"{R26}/tables/per_kinase_support.tsv", sep="\t", float_format="%.3f")

# strata
rows = []
for var in ["evidence", "LT_bin", "motif_bin", "mgroup"]:
    for lv, d in U.groupby(var, observed=True):
        rows.append(dict(variable=var, level=lv, n=len(d), support_rate=d.supported.mean(),
                         unsupported_rate=(d.cls == "unsupported").mean(), median_mean_z=d.mean_z.median()))
S = pd.DataFrame(rows); S.to_csv(f"{R26}/tables/relations_strata.tsv", sep="\t", index=False, float_format="%.3f")
print(S.round(3).to_string())

# models (group fixed effects; groups with >=10 testable units)
big = pg.index[pg.n_testable >= 10]
M = U[U.group.isin(big)].copy()
out = []
for f in ["supported ~ C(group) + in_vivo + in_vitro", "supported ~ C(group) + logLT + logMS",
          "supported ~ C(group) + motif", "supported ~ C(group) + in_vivo + in_vitro + logLT + logMS + motif"]:
    d = M.dropna(subset=["motif"]) if "motif" in f else M
    r = smf.logit(f, d).fit(disp=0)
    for v in r.params.index:
        if not v.startswith("C(") and v != "Intercept":
            out.append(dict(model=f, term=v, OR=np.exp(r.params[v]), lo=np.exp(r.conf_int().loc[v, 0]),
                            hi=np.exp(r.conf_int().loc[v, 1]), p=r.pvalues[v], n=int(r.nobs)))
for f in ["mean_z ~ C(group) + in_vivo + in_vitro + logLT + logMS + motif"]:
    r = smf.ols(f, M.dropna(subset=["motif"])).fit()
    for v in r.params.index:
        if not v.startswith("C(") and v != "Intercept":
            out.append(dict(model=f, term=v, OR=r.params[v], lo=r.conf_int().loc[v, 0], hi=r.conf_int().loc[v, 1], p=r.pvalues[v], n=int(r.nobs)))
R = pd.DataFrame(out); R.to_csv(f"{R26}/tables/relations_models.tsv", sep="\t", index=False, float_format="%.4g")
print(R.round(3).to_string())
for v in ["motif", "LT", "MS"]:
    d = M.dropna(subset=[v])
    print(v, "spearman with mean_z (within-group pooled ranks):",
          spearmanr(d.groupby("group")[v].rank(pct=True), d.groupby("group").mean_z.rank(pct=True)))
