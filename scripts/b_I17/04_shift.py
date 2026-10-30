"""Potency shift dpEC50 = pEC50_site - pKd_target (Klaeger Kd_app) and its determinants.
Target set T(drug) = Klaeger kinases with Kd <= 3x the most potent designated target (on-target cluster).
Annotations per (experiment, site): self (site on a T protein), PSP/OmniPath substrate of T, predicted substrate
(kinaseSub top5 contains a T kinase), STRING(>=700) distance from T, COMPARTMENTS localisation, direction.
Outputs: results/b_I17/shift_sites.tsv, shift_by_experiment.tsv, shift_models.tsv
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import re, numpy as np, pandas as pd, networkx as nx
import statsmodels.formula.api as smf
OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"

F = pd.read_parquet(f"{RES}/curves_all.parquet")
E = pd.read_csv(f"{OUT}/experiments.tsv", sep="\t").set_index("exp_id")
K = pd.read_parquet(f"{OUT}/klaeger_kd.parquet")
DM = pd.read_csv(f"{RES}/drug_map.tsv", sep="\t").set_index("agent")
s1 = pd.read_excel("/workspace/data/Klaeger/aan4368_Table_S1.xlsx", sheet_name="Inhibitor annotation").set_index("Drug")
S = pd.read_parquet(f"{OUT}/sites_min.parquet").set_index("site")
KS = pd.read_parquet(f"{OUT}/ks_edges.parquet")
P5 = pd.read_parquet(f"{OUT}/pred_top5.parquet")
ST = pd.read_parquet(f"{OUT}/string_edges.parquet").dropna()
loc = pd.read_csv(f"{OUT}/loc.tsv", sep="\t", index_col=0)

ks_by_site = KS.groupby("site").kinase.apply(set)
pred_by_site = P5[P5["rank"] <= 5].groupby("site_id").gene.apply(set)
G = nx.Graph(); G.add_edges_from(zip(ST.a, ST.b))


def drug_targets(agent):
    if agent not in DM.index or not DM.loc[agent, "all_mapped"] or ";" in agent:
        return None
    kd = DM.loc[agent, "klaeger"]
    des = [t.strip() for t in str(s1.loc[kd, "Designated targets"]).split(",")]
    col = K[kd]
    dk = col.reindex(des).dropna(); dk = dk[np.isfinite(dk)]
    if len(dk) == 0:
        return None
    kd0 = dk.min(); tgt = dk.idxmin()
    fin = col[np.isfinite(col)].dropna()
    T = set(g for x in fin[fin <= 3 * kd0].index for g in str(x).split(";"))
    return dict(klaeger=kd, target=tgt, kd_target=kd0, pKd_target=9 - np.log10(kd0), T=T, profile=fin)


rows, exrows = [], []
for ex, g in F[F.sig].groupby("exp_id"):
    e = E.loc[ex]
    if e.pclass != "kinase_inhibitor": continue
    dt = drug_targets(e.agent)
    if dt is None: continue
    T = dt["T"]
    dist = {}
    for t in T:
        if t in G:
            for n, d in nx.single_source_shortest_path_length(G, t, cutoff=4).items():
                dist[n] = min(d, dist.get(n, 9))
    g = g.copy()
    g["gene"] = S.gene.reindex(g.site).values
    g["self"] = g.gene.isin(T)
    g["ks_direct"] = [bool(ks_by_site.get(s, set()) & T) for s in g.site]
    g["ks_any"] = g.site.isin(ks_by_site.index)
    g["pred_direct"] = [bool(pred_by_site.get(s, set()) & T) for s in g.site]
    g["string_dist"] = [0 if gn in T else dist.get(gn, 5) if isinstance(gn, str) else np.nan for gn in g.gene]
    g["nucleus"] = loc.reindex(g.gene)["Nucleus"].fillna(0).values >= 4
    g["pmembrane"] = loc.reindex(g.gene)["Plasma membrane"].fillna(0).values >= 4
    g["direction"] = np.where(g.delta < 0, "down", "up")
    g["pKd_target"] = dt["pKd_target"]; g["target"] = dt["target"]; g["klaeger"] = dt["klaeger"]
    g["shift"] = g.pEC50 - dt["pKd_target"]
    # dose-window: maximal shift observable
    doses = np.array([float(v) for v in e.doses.split(";")]); nz = doses[doses > 0]
    g["shift_min_obs"] = 9 - np.log10(nz.max()) - dt["pKd_target"]
    g["shift_max_obs"] = 9 - np.log10(nz.min()) - dt["pKd_target"]
    rows.append(g)
    d = g[g.direction == "down"]
    exrows.append(dict(exp_id=ex, agent=e.agent, cell_line=e.cell_line, time_min=e.time_min, target=dt["target"],
                       kd_target_nM=dt["kd_target"], pKd_target=dt["pKd_target"], n_T=len(T), n_sig=len(g), n_down=len(d),
                       n_up=(g.direction == "up").sum(), median_pEC50_down=d.pEC50.median(),
                       median_shift_down=d["shift"].median(), iqr_shift_down=d["shift"].quantile(.75) - d["shift"].quantile(.25),
                       frac_within_0p5=(d["shift"].abs() <= 0.5).mean(), frac_amplified=(d["shift"] > 1).mean(),
                       frac_attenuated=(d["shift"] < -1).mean(), n_direct_down=(d.ks_direct | d.self).sum(),
                       median_shift_direct=d.loc[d.ks_direct | d.self, "shift"].median(),
                       median_hill_down=d.hill.median(), window_lo=g.shift_min_obs.iloc[0], window_hi=g.shift_max_obs.iloc[0]))
A = pd.concat(rows, ignore_index=True)
X = pd.DataFrame(exrows)
A.to_csv(f"{RES}/shift_sites.tsv", sep="\t", index=False)
X.to_csv(f"{RES}/shift_by_experiment.tsv", sep="\t", index=False)
print(X.round(2).to_string())

# models: experiment fixed effects, cluster-robust SE by drug
A["log2hill"] = np.log2(A.hill)
A["direct"] = A.ks_direct | A.self
A["dist_c"] = A.string_dist.clip(upper=4)
A["down"] = (A.direction == "down").astype(int)
A["drug"] = A.agent
mres = []
for y in ["shift", "log2hill"]:
    for form, sub in [("direct + pred_direct + C(dist_c) + nucleus + pmembrane + down", A),
                      ("direct + pred_direct + C(dist_c) + nucleus + pmembrane", A[A.down == 1])]:
        d = sub.dropna(subset=["dist_c"])
        m = smf.ols(f"{y} ~ {form} + C(exp_id)", data=d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.drug)[0]})
        for term in m.params.index:
            if term.startswith("C(exp_id)") or term == "Intercept": continue
            mres.append(dict(outcome=y, subset="all" if len(sub) == len(A) else "down_only", term=term, coef=m.params[term],
                             se=m.bse[term], p=m.pvalues[term], n=int(m.nobs), n_drugs=d.drug.nunique()))
M = pd.DataFrame(mres); M.to_csv(f"{RES}/shift_models.tsv", sep="\t", index=False)
print(M.round(4).to_string())
