"""Task 1b: can AKT vs RSK be told apart? Motif percentiles; per-site AKT-verified conditions where site did not fall;
AKT-marker behaviour in MAPK-verified conditions; RSK-genetic (Roux RNAi) values."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__)); from lib_f import *
OUT = "/workspace/results/audit_F/tables"
pd.set_option("display.width", 250)
SITES = ["O14757_S280", "P29317_S897", "Q14247_S405", "P36507_T394", "P26358_S143", "O00443_S259", "Q01167_S373", "Q96F86_S161",
         "Q9UJU6_S269", "P0C1Z6_S180", "Q9UKA4_T1100"]
mot = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
KS = ["AKT1", "AKT2", "AKT3", "SGK1", "RPS6KB1", "RPS6KA1", "RPS6KA3", "RPS6KA5", "PIM1", "PRKACA", "PRKCA", "PRKCD", "MAPK1", "MAPK3",
      "CDK1", "CDK2", "CDK5", "MTOR", "GSK3B"]
M = mot.reindex(SITES)[[k for k in KS if k in mot.columns]].astype(float)
SA = load.sites().set_index("site")
M.insert(0, "window", SA.sequence_window.reindex(SITES)); M.insert(0, "gene", SA.gene.reindex(SITES))
M["AKT_max"] = M[["AKT1", "AKT2", "AKT3"]].max(1); M["RSK_max"] = M[["RPS6KA1", "RPS6KA3"]].max(1)
M["RSK_minus_AKT"] = M.RSK_max - M.AKT_max
M.to_csv(f"{OUT}/motif_percentiles.tsv", sep="\t", float_format="%.0f")
print(M.round(0).to_string())
# genome-wide: how correlated are AKT and RSK percentiles (can motif separate them?)
mm = mot[["AKT1", "RPS6KA1", "RPS6KA3", "MAPK1"]].dropna()
print("Spearman AKT1 vs RPS6KA1 all sites:", round(mm.AKT1.corr(mm.RPS6KA1, method="spearman"), 3),
      " vs RPS6KA3:", round(mm.AKT1.corr(mm.RPS6KA3, method="spearman"), 3))
hi = mm[(mm.RPS6KA1 >= 90) | (mm.AKT1 >= 90)]
print("Among sites >=90 for AKT1 or RSK1: fraction >=90 for both:", round(((hi.RPS6KA1 >= 90) & (hi.AKT1 >= 90)).mean(), 3), "n", len(hi))

L = pd.read_csv(f"{OUT}/site_condition_long.tsv", sep="\t", dtype={"study": str})
rows = []
for s in SITES:
    d = L[L.site == s]
    a = d[(d.cls == "AKT") & d.verified]
    m = d[d.cls.isin(["MEK", "ERK", "RAF"]) & d.verified & ~d.lab.eq("Thibault")]
    rows.append(dict(gene=M.gene[s], site=s, AKTi_verified_conds=len(a), AKTi_labs=";".join(sorted(a.lab.unique())),
                     AKTi_conds_site_not_falling=int((a.z > -1).sum()), AKTi_median_site_z=a.z.median(), AKTi_median_AKTmark=a.marker_z.median(),
                     AKTi_median_MAPKmark=a.MAPK_mark.median(),
                     MAPKi_verified_conds=len(m), MAPKi_median_site_z=m.z.median(), MAPKi_median_AKTmark=m.AKT_mark.median(),
                     MAPKi_conds_AKTmark_le_m2=int((m.AKT_mark <= -2).sum())))
B = pd.DataFrame(rows)
B.to_csv(f"{OUT}/akt_vs_mapk_verification.tsv", sep="\t", index=False, float_format="%.2f")
print(B.round(2).to_string())
# detailed AKT-verified condition list
det = L[(L.cls == "AKT") & L.verified][["gene", "conditionID", "lab", "agent", "cell", "time_min", "z", "marker_z", "MAPK_mark"]]
det.to_csv(f"{OUT}/akt_verified_conditions_per_site.tsv", sep="\t", index=False, float_format="%.2f")
# Roux RSK genetic + BI-D1870
C = load.conditions()
rc = C[C.study_pmid.eq("25002506")].conditionID.tolist()
X = load.sites_matrix(rc)
med = X.median(); mad = (X - med).abs().median() * 1.4826
Zr = ((X - med) / mad.clip(lower=0.05)).clip(-8, 8)
lab = dict(zip(C.conditionID, C.condition_raw + " | " + C.cell_line_norm.astype(str)))
R = Zr.reindex(SITES + ["Q07889_S1161", "P49841_S9", "Q96B36_T246"]).rename(columns=lab).round(2)
R.insert(0, "gene", SA.gene.reindex(R.index))
R.to_csv(f"{OUT}/roux_rsk_conditions.tsv", sep="\t")
print(R.to_string())
