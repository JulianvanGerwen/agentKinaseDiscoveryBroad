"""I26 step 4: for powered-unsupported PSP edges, find alternative inhibitor groups whose LoF consistently reduces the site.
Alternative H must be 'supported' for the site (>=2 verified-inhibition studies with study-mean z<=-1, >=50%).
Also run with no positive-control filter (L_none) to let groups without verified inhibition (e.g. CK2) contribute (flagged).
Output tables/corrections.tsv (one row per unsupported edge x alternative), tables/corrections_top.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *

E, Z, GC = load_base()
gold = gold_sets(E)
Lm = pd.read_parquet(f"{C26}/L_main.parquet"); Lm["group"] = Lm.group.astype(str); Lm["study"] = Lm.study.astype(str)
Ln, _ = study_long(Z, GC, gold, thr=None)
Ln.to_parquet(f"{C26}/L_none.parquet")
U = pd.read_csv(f"{C26}/units_testable.tsv", sep="\t")
uns = U[U.cls == "unsupported"]
sites = set(uns.site)
Am = classify(Lm[Lm.site.isin(sites)]); An = classify(Ln[Ln.site.isin(sites)])
G = pd.read_csv(f"{C26}/groups.tsv", sep="\t"); members = G.groupby("group").kinase.apply(list).to_dict()
mot = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
DOWN = {"MEK": ["MAPK1", "MAPK3", "RPS6KA1", "RPS6KA3", "RPS6KA5", "MKNK1"],
        "RAF": ["MAPK1", "MAPK3", "RPS6KA1", "RPS6KA3", "RPS6KA5", "MKNK1"],
        "EGFR/ERBB2": ["MAPK1", "MAPK3", "RPS6KA1", "RPS6KA3", "AKT1", "RPS6KB1"],
        "MET": ["MAPK1", "MAPK3", "RPS6KA1", "AKT1"], "ALK/ROS1": ["MAPK1", "MAPK3", "RPS6KA1", "AKT1"],
        "ABL/KIT/PDGFR": ["MAPK1", "MAPK3", "AKT1"], "MTOR": ["RPS6KB1", "RPS6KB2", "SGK1"], "AKT": ["SGK1"],
        "ERK": ["RPS6KA1", "RPS6KA3", "RPS6KA5", "MKNK1"]}
AXIS = {"EGFR/ERBB2": "MAPK", "MET": "MAPK", "ALK/ROS1": "MAPK", "ABL/KIT/PDGFR": "MAPK", "RAF": "MAPK", "MEK": "MAPK", "ERK": "MAPK",
        "AKT": "PI3K", "MTOR": "PI3K", "S6K": "PI3K"}
UPSTREAM = {"EGFR/ERBB2", "MET", "ALK/ROS1", "ABL/KIT/PDGFR", "RAF", "MEK"}
KS = pd.read_csv("/workspace/results/infra_compendium/psp_ks_edges_human_mapped.tsv", sep="\t")
psp_by_site = KS.groupby("site").kinase.apply(lambda x: ";".join(sorted(set(x))))
SA = pd.read_parquet("/workspace/cache/b_I03/site_annot.parquet").set_index("site")


def best_motif(site, ks):
    vals = [(mot.at[site, k], k) for k in ks if k in mot.columns and site in mot.index]
    return max(vals) if vals else (np.nan, None)


rows = []
for r in uns.itertuples():
    s, g0 = r.site, r.group
    for tag, A in [("verified", Am), ("unverified_ctrl", An)]:
        a = A[(A.site == s) & (A.group != g0) & (A.cls == "supported")]
        n_groups_sup = ((A.site == s) & (A.cls == "supported")).sum()
        for h in a.itertuples():
            mp, mk = best_motif(s, members[h.group])
            dp, dk = best_motif(s, DOWN.get(h.group, []))
            mo, _ = best_motif(s, members[g0])
            others = A[(A.site == s) & (A.group != h.group) & (A.n_meas >= 2)]
            rows.append(dict(site=s, gene=r.gene, psp_group=g0, psp_kinases=r.kinases, psp_mean_z=r.mean_z, psp_n_meas=r.n_meas,
                             psp_motif=mo, alt_group=h.group, evidence=tag, alt_n_meas=h.n_meas, alt_n_dec=h.n_dec, alt_mean_z=h.mean_z,
                             alt_motif_member=mk, alt_motif_pct=mp, downstream_motif_kinase=dk, downstream_motif_pct=dp,
                             n_groups_supporting_site=n_groups_sup, other_groups_median_z=others.mean_z.median(),
                             same_axis=AXIS.get(h.group, h.group) == AXIS.get(g0, g0), alt_is_upstream=h.group in UPSTREAM,
                             all_psp_kinases_site=psp_by_site.get(s, ""), LT_LIT=r.LT, MS_LIT=r.MS, evidence_psp=r.evidence,
                             window=SA.sequence_window.get(s, "")))
R = pd.DataFrame(rows)
R["alt_already_psp"] = [any(k in str(p).split(";") for k in members[h]) for h, p in zip(R.alt_group, R.all_psp_kinases_site)]
R["specificity"] = R.alt_mean_z - R.other_groups_median_z
R["score"] = (-R.alt_mean_z).clip(upper=6) * np.minimum(R.alt_n_dec, 4) / 4 + 0.02 * R[["alt_motif_pct", "downstream_motif_pct"]].max(axis=1).fillna(50) \
    - 0.5 * (R.n_groups_supporting_site - 1).clip(lower=0) - 1.0 * (R.evidence != "verified")
R = R.sort_values("score", ascending=False)
R.to_csv(f"{R26}/tables/corrections.tsv", sep="\t", index=False, float_format="%.3g")
print(len(uns), "unsupported units;", R.site.nunique(), "with >=1 alternative;",
      R[R.evidence == "verified"].site.nunique(), "verified-only")
print(R.groupby(["psp_group", "alt_group", "evidence"]).size().sort_values(ascending=False).head(30).to_string())
# top list: one row per (site, psp_group), best alternative; exclude same-axis upstream (pathway-consistent) for misattribution list
best = R.drop_duplicates(["site", "psp_group"])
best.to_csv(f"{R26}/tables/corrections_best_per_edge.tsv", sep="\t", index=False, float_format="%.3g")
cols = ["gene", "site", "psp_kinases", "psp_mean_z", "psp_n_meas", "psp_motif", "alt_group", "evidence", "alt_n_dec", "alt_n_meas",
        "alt_mean_z", "alt_motif_member", "alt_motif_pct", "downstream_motif_kinase", "downstream_motif_pct", "n_groups_supporting_site",
        "same_axis", "alt_already_psp", "LT_LIT", "MS_LIT", "evidence_psp", "window"]
print(best[cols].head(40).round(2).to_string())
