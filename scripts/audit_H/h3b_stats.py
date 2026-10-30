"""Audit H3 stats: fast vs stable coupling, adjusted (sd, missingness, n_tum), matched, site-cluster bootstrap,
own-minus-other-family specificity, non-gold / no-same-protein-gold restriction."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, statsmodels.formula.api as smf
sys.path.insert(0, "/workspace/scripts/b_I26"); from lib26 import gold_sets
OUT = "/workspace/results/audit_H"; rng = np.random.default_rng(9)
D = pd.read_csv(f"{OUT}/h3_site_cohort.tsv", sep="\t")
E = pd.read_parquet("/workspace/cache/b_I26/psp_edges_annot.parquet"); gold = gold_sets(E)
gprot = {g: {s.split("_")[0] for s in v} for g, v in gold.items()}
D["prot"] = D.site.str.split("_").str[0]
D["prot_has_gold"] = [p in gprot.get(g, ()) for g, p in zip(D.group, D.prot)]
D["spec"] = D.coupling - D.coupling_other
S = D[(D.kind == "substrate") & D.cls.isin(["fast", "stable"])].dropna(subset=["coupling"]).copy()
S["fast"] = (S.cls == "fast").astype(float); S["lognt"] = np.log(S.n_tum)
rows = []


def boot_diff(d, y, B=500):
    """site-cluster bootstrap of adjusted fast-stable difference: OLS y ~ fast + sd + miss + lognt + C(cohort:group)."""
    f = f"{y} ~ fast + sd_resid + miss + lognt + C(cg)"
    d = d.assign(cg=d.cohort + d.group).dropna(subset=[y])
    est = smf.ols(f, d).fit().params.fast
    sites = d.site.unique(); pos = d.groupby("site").indices
    bs = []
    for _ in range(B):
        ix = np.concatenate([pos[s] for s in rng.choice(sites, len(sites))])
        bs.append(smf.ols(f, d.iloc[ix]).fit().params.fast)
    return est, np.percentile(bs, 2.5), np.percentile(bs, 97.5), d.site.nunique(), len(d)


def matched(d, y):
    """within cohort x group: each fast matched to nearest stable on (sd_resid, miss) ranks; mean diff; site-cluster bootstrap."""
    pairs = []
    for _, x in d.groupby(["cohort", "group"]):
        f = x[x.fast == 1]; s = x[x.fast == 0]
        if len(f) < 3 or len(s) < 3: continue
        Z = x[["sd_resid", "miss"]].rank(pct=True)
        sf, ss = Z.loc[f.index].values, Z.loc[s.index].values
        dist = ((sf[:, None, :] - ss[None, :, :]) ** 2).sum(-1)
        j = dist.argmin(1)
        pairs.append(pd.DataFrame({"site": f.site.values, "diff": f[y].values - s[y].values[j],
                                   "caliper": np.sqrt(dist.min(1))}))
    P = pd.concat(pairs); P = P[P.caliper < 0.1].dropna()
    sites = P.site.unique(); pos = P.groupby("site").indices
    bs = [P["diff"].values[np.concatenate([pos[s] for s in rng.choice(sites, len(sites))])].mean() for _ in range(1000)]
    return P["diff"].mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(sites), len(P)


for v in ["I28", "isch", "isch+5PC"]:
    for sub, mask in [("all", None), ("non-gold", "~gold"), ("no gold on same protein", "~prot_has_gold")]:
        d = S[S.variant == v]
        if mask == "~gold": d = d[~d.gold.astype(bool)]
        if mask == "~prot_has_gold": d = d[~d.prot_has_gold]
        for y in ["coupling", "spec"]:
            med_f = d[d.fast == 1].groupby("site")[y].mean().median(); med_s = d[d.fast == 0].groupby("site")[y].mean().median()
            e, lo, hi, ns, n = boot_diff(d, y, B=200)
            me, mlo, mhi, mns, mn = matched(d, y)
            rows.append(dict(variant=v, subset=sub, outcome=y, site_median_fast=med_f, site_median_stable=med_s,
                             adj_diff=e, adj_lo=lo, adj_hi=hi, n_sites=ns, n_obs=n, matched_diff=me, m_lo=mlo, m_hi=mhi, m_sites=mns))
        print(v, sub, flush=True)
O = pd.DataFrame(rows); O.to_csv(f"{OUT}/h3_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(O.round(3).to_string())
