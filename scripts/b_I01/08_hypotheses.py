"""Candidate hypotheses for dark kinases with reliable kinase-specific signal.
Reliability FDR from normal approx. of rel_z (50 group-matched random-set nulls, protein-split halves).
For each reliable dark kinase: top non-family partners by a_resid (row z-scored), with orthogonal support
(DepMap r, GTEx rho, PPI exp/lit, HPA Jaccard, COMPARTMENTS cosine) and top perturbation contexts."""
import numpy as np, pandas as pd, sys, os
from pathlib import Path
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
OUT = Path("/workspace/results/b_I01"); CA = OUT / "cache"
T = pd.read_csv(OUT / "tables/specificity_per_kinase.tsv", sep="\t", index_col=0)
M = pd.read_csv(OUT / "tables/substrate_sets.tsv", sep="\t", index_col=0)
kin = list(M.index)
T["rel_pz"] = norm.sf(T.rel_z)
for g in [True, False]:
    m = (T.dark == g) & T.rel_pz.notna()
    T.loc[m, "rel_fdr"] = multipletests(T.loc[m, "rel_pz"], method="fdr_bh")[1]
T.to_csv(OUT / "tables/specificity_per_kinase.tsv", sep="\t")
print("reliable (FDR<=0.05): dark", ((T.rel_fdr <= 0.05) & T.dark).sum(), "/", (T.dark & T.rel_pz.notna()).sum(),
      "; non-dark", ((T.rel_fdr <= 0.05) & ~T.dark).sum(), "/", (~T.dark & T.rel_pz.notna()).sum())
A = pd.read_parquet(CA / "assoc_resid.parquet").loc[kin, kin]
Af = pd.read_parquet(CA / "assoc_full.parquet").loc[kin, kin]
Ac = pd.read_parquet(CA / "assoc_conf.parquet").loc[kin, kin]
L = {n: pd.read_parquet(CA / f"lab_{n}.parquet").reindex(index=kin, columns=kin) for n in ["depmap", "gtex", "ppi_exp", "ppi_lit", "hpa", "compartments"]}
dep_thr = np.nanquantile(L["depmap"].values[np.triu_indices(len(kin), 1)], 0.98)
W = pd.read_csv(OUT / "tables/contexts_wide.tsv", sep="\t")
rel = T[(T.rel_fdr <= 0.05) & T.dark].index
rows = []
for k in rel:
    excl = set(M.loc[k, "neighbours"].split(";")) | set(M.index[M.family == M.loc[k, "family"]])
    a = A.loc[k].drop(list(excl), errors="ignore").dropna()
    z = (a - a.mean()) / a.std()
    for j in z.sort_values(ascending=False).index[:10]:
        rows.append(dict(dark_kinase=k, family=M.loc[k, "family"], partner=j, partner_family=M.loc[j, "family"],
                         a_resid=A.loc[k, j], a_resid_z=z[j], a_full=Af.loc[k, j], a_conf=Ac.loc[k, j],
                         depmap_r=L["depmap"].loc[k, j], depmap_top2pct=L["depmap"].loc[k, j] >= dep_thr,
                         gtex_rho=L["gtex"].loc[k, j], ppi_exp=L["ppi_exp"].loc[k, j], ppi_lit=L["ppi_lit"].loc[k, j],
                         hpa_jacc=L["hpa"].loc[k, j], comp_cos=L["compartments"].loc[k, j]))
P = pd.DataFrame(rows)
P["n_orth_support"] = P.depmap_top2pct.astype(int) + (P.ppi_exp > 0).astype(int) + (P.ppi_lit > 0).astype(int) + (P.gtex_rho > 0.7).astype(int)
P.to_csv(OUT / "tables/dark_partners_top10.tsv", sep="\t", index=False)
# depmap co-dependency partners of reliable dark kinases: what is the a_resid rank of the top-codependent kinase?
cd = []
for k in rel:
    excl = set(M.loc[k, "neighbours"].split(";")) | set(M.index[M.family == M.loc[k, "family"]])
    d = L["depmap"].loc[k].drop(list(excl), errors="ignore").dropna()
    a = A.loc[k].reindex(d.index)
    for j in d[d >= dep_thr].index:
        if pd.isna(a[j]): continue
        cd.append(dict(dark_kinase=k, partner=j, depmap_r=d[j], a_resid=a[j], a_resid_pct=(a.dropna() < a[j]).mean(),
                       a_full_pct=(Af.loc[k].reindex(d.index).dropna() < Af.loc[k, j]).mean(),
                       a_conf_pct=(Ac.loc[k].reindex(d.index).dropna() < Ac.loc[k, j]).mean()))
CD = pd.DataFrame(cd).sort_values("a_resid_pct", ascending=False)
CD.to_csv(OUT / "tables/dark_codependency_partners.tsv", sep="\t", index=False)
Cx = W[W.kinase.isin(rel) & (W.p_resid <= 0.01) & (W.frac_same_sign >= 0.8) & (W.n_study >= 3)].sort_values("p_resid")
Cx.to_csv(OUT / "tables/dark_contexts_top.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(T.loc[rel, ["family", "n_spec", "r2_full_conf", "r2_spec_conf", "rel_obs", "rel_null_mean", "rel_z", "rel_fdr"]].sort_values("rel_z", ascending=False).round(3).to_string())
print(P[P.n_orth_support >= 1].round(3).to_string())
print(CD.head(30).round(3).to_string())
print(Cx.round(4).to_string())
