"""I28 step 6: are fast-turnover sites 'reporter' sites in tumours?
Per CPTAC cohort (tumours only): site log2 residualised on host-protein abundance (OLS per site, >=20 tumours).
(1) variability: SD of raw and protein-residualised values; (2) coupling: Spearman of residual with the family's
leave-one-out footprint = mean of within-cohort z-scored residuals of the family's other gold (PSP in-vivo) substrates.
Fast vs stable compared within family, per cohort, adjusted for n tumours measured & mean intensity (OLS), meta over cohorts.
Matched random sites provide the baseline coupling.
Out: tables/cptac_site_stats.tsv, tables/cptac_fast_vs_stable.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
import statsmodels.api as sm
sys.path.insert(0, os.path.dirname(__file__)); from lib28 import *
sys.path.insert(0, "/workspace/scripts/b_I26"); from lib26 import gold_sets
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
A = pd.read_csv(f"{R28}/tables/turnover_atlas.tsv", sep="\t")
L = pd.read_parquet(f"{C28}/L_units.parquet")
from lib28 import classify_units
AR = classify_units(L[L.set == "random"]).assign(source="random")
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]
U = pd.concat([A[A.group.isin(MAIN) & (A.n_units >= 2)][["group", "site", "cls", "R"]].assign(kind="substrate"),
               AR[AR.group.isin(MAIN) & (AR.n_units >= 2)][["group", "site", "cls", "R"]].assign(kind="random")])
E = pd.read_parquet(f"{C26}/psp_edges_annot.parquet"); gold = gold_sets(E)
rows = []
for co in cl.COHORTS:
    t = cl.tumour_ids(co)
    ph = cl.matrix(co, "phospho"); ph = ph[[c for c in ph.columns if c in t]]
    pr = cl.matrix(co, "protein"); pr = pr.reindex(columns=ph.columns)
    need = set(U.site) | set().union(*[gold[g] for g in MAIN])
    ph = ph[ph.index.isin(need)]
    ph = ph[ph.notna().sum(axis=1) >= 20]
    res = {}
    for s, v in ph.iterrows():
        p = pr.loc[s.split("_")[0]] if s.split("_")[0] in pr.index else None
        y = v.values.astype(float)
        if p is not None:
            x = p.values.astype(float); ok = ~np.isnan(y) & ~np.isnan(x)
            if ok.sum() >= 20:
                b = np.polyfit(x[ok], y[ok], 1); r = np.full_like(y, np.nan); r[ok] = y[ok] - np.polyval(b, x[ok]); res[s] = r; continue
        res[s] = np.full_like(y, np.nan)
    R = pd.DataFrame(res, index=ph.columns).T
    Rz = R.sub(R.mean(axis=1), axis=0).div(R.std(axis=1), axis=0)
    for g in MAIN:
        gs = [s for s in gold[g] if s in Rz.index and Rz.loc[s].notna().sum() >= 20]
        if len(gs) < 4: continue
        Gs = Rz.loc[gs]; tot = Gs.sum(); cnt = Gs.notna().sum()
        for _, u in U[U.group == g].iterrows():
            s = u.site
            if s not in ph.index: continue
            if s in gs:
                fp = (tot - Gs.loc[s].fillna(0)) / (cnt - Gs.loc[s].notna())
            else:
                fp = tot / cnt
            r = R.loc[s]; ok = r.notna() & fp.notna()
            coup = r[ok].rank().corr(fp[ok].rank()) if ok.sum() >= 20 else np.nan
            rows.append(dict(cohort=co, group=g, site=s, kind=u.kind, cls=u.cls, R=u.R, n_tum=int(ph.loc[s].notna().sum()),
                             mean_raw=ph.loc[s].mean(), sd_raw=ph.loc[s].std(), sd_resid=r.std(), coupling=coup, n_gold=len(gs)))
    print(co, len(rows), flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{R28}/tables/cptac_site_stats.tsv", sep="\t", index=False)
out = []
for outc in ["sd_raw", "sd_resid", "coupling"]:
    for comp in [("fast", "stable", "substrate"), ("fast", "random_all", None)]:
        ests = []
        for (co, g), d in D.groupby(["cohort", "group"]):
            if comp[1] == "stable":
                d = d[(d.kind == "substrate") & d.cls.isin(["fast", "stable"])]; x = (d.cls == "fast").astype(float)
            else:
                d = d[((d.kind == "substrate") & (d.cls == "fast")) | (d.kind == "random")]; x = (d.kind == "substrate").astype(float)
            d = d.dropna(subset=[outc, "mean_raw"])
            if x.loc[d.index].sum() < 5 or (1 - x.loc[d.index]).sum() < 5: continue
            Xm = sm.add_constant(pd.DataFrame({"x": x.loc[d.index], "lognt": np.log(d.n_tum), "mean": d.mean_raw}))
            m = sm.OLS(d[outc], Xm).fit()
            ests.append((m.params.x, m.bse.x, co, g, len(d)))
        if not ests: continue
        e = np.array([a[0] for a in ests]); se = np.array([a[1] for a in ests]); w = 1 / se ** 2
        b = (w * e).sum() / w.sum(); sb = 1 / np.sqrt(w.sum())
        from scipy.stats import norm
        out.append(dict(outcome=outc, contrast=f"{comp[0]}_vs_{comp[1]}", n_strata=len(ests), beta=b, lo=b - 1.96 * sb, hi=b + 1.96 * sb,
                        p=2 * norm.sf(abs(b / sb)), frac_strata_pos=(e > 0).mean()))
O = pd.DataFrame(out); O.to_csv(f"{R28}/tables/cptac_fast_vs_stable.tsv", sep="\t", index=False)
print(O.round(4).to_string())
print(D.groupby(["kind", "cls"]).coupling.median().round(3))
# continuous: coupling vs -R within group (substrates)
from scipy.stats import spearmanr
x = D[D.kind == "substrate"].groupby(["group", "site"]).agg(coupling=("coupling", "mean"), R=("R", "first")).dropna()
print("site-level mean coupling vs -R, within-group spearman:", {g: round(spearmanr(d.coupling, -d.R)[0], 3) for g, d in x.groupby(level=0) if len(d) > 20})
