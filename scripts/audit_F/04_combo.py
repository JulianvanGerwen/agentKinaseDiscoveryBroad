"""Task 1b extra: decryptM A549 30 min dose series: MK2206 vs selumetinib vs selumetinib+MK2206 (same lab, same cell, same time).
If a site is a dual AKT/RSK site masked by feedback, the combination should lower it more than selumetinib alone."""
import sys; sys.path.insert(0, "/workspace/scripts/infra_compendium")
import numpy as np, pandas as pd, load
pd.set_option("display.width", 250)
C = load.conditions()
c = C[(C.study_pmid == "36926954") & (C.cell_line_norm == "A549") & C.dose_series_id.isin(["DS0011", "DS0019", "DS0010"])]
c = c[c.dose_nM >= 300]  # top doses (single >=300 nM; combo total >=300 nM)
X = load.sites_matrix(c.conditionID.tolist())
SITES = {"O14757_S280": "CHEK1 S280", "P29317_S897": "EPHA2 S897", "Q14247_S405": "CTTN S405", "P36507_T394": "MAP2K2 T394",
         "P26358_S143": "DNMT1 S143", "O00443_S259": "PIK3C2A S259", "Q01167_S373": "FOXK2 S373", "Q96F86_S161": "EDC3 S161",
         "Q9UJU6_S269": "DBNL S269", "P0C1Z6_S180": "TFPT S180", "Q9UKA4_T1100": "AKAP11 T1100",
         "P49841_S9": "[AKT/RSK] GSK3B S9", "Q96B36_T246": "[AKT] PRAS40 T246", "P49815_T1462": "[AKT] TSC2 T1462", "O43524_S253": "[AKT] FOXO3 S253",
         "P28482_Y187": "[MEK] ERK2 Y187", "P27361_Y204": "[MEK] ERK1 Y204", "Q07889_S1161": "[RSK] SOS1 S1161", "P28482_T185": "[MEK] ERK2 T185"}
lab = {"DS0011": "MK2206", "DS0019": "selumetinib", "DS0010": "selu+MK2206"}
out = {}
for ds, d in c.groupby("dose_series_id"):
    sub = X.reindex(list(SITES))[d.conditionID]
    out[lab[ds] + " mean log2FC"] = sub.mean(axis=1); out[lab[ds] + " n"] = sub.notna().sum(axis=1)
R = pd.DataFrame(out); R.index = [SITES[s] for s in R.index]
R["combo_minus_selu"] = R["selu+MK2206 mean log2FC"] - R["selumetinib mean log2FC"]
R.round(2).to_csv("/workspace/results/audit_F/tables/decryptM_A549_combo.tsv", sep="\t")
print(R.round(2).to_string())
