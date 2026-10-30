"""Map compendium drug x cell to GDSC1/GDSC2/PRISM AUC + ln IC50. Out: tables/T2_sens_matched.tsv, T0_drug_map.tsv"""
import re, pandas as pd, numpy as np
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())
L = pd.read_parquet("/workspace/cache/b_I23/sens_long.parquet")
F = pd.read_csv("/workspace/results/b_I23/tables/T1_features.tsv", sep="\t")
# vetted synonyms (compendium name -> list of screen names); AZD1480's GDSC synonym hit on Selumetinib(1736) is an annotation error -> not used
SYN = {"JNKIN8": ["JNK Inhibitor VIII"], "ceralasertib": ["AZD6738"], "capivasertib": ["AZD5363"], "ipatasertib": ["GDC-0068"],
       "ulixertinib": ["BVD-523"], "belumosudil": ["KD025"], "torkinib": ["PP242"], "torin1": ["torin-1"], "TAK715": ["TAK-715"],
       "MK2206": ["MK-2206"], "PF4708671": ["PF-4708671"], "VE821": ["VE-821"], "pictilisib": ["GDC-0941"], "KU0063794": ["ku-0063794"],
       "BIRB796": ["Doramapimod"], "TAE684": ["NVP-TAE684"], "CX4945": ["silmitasertib"], "PD184352": ["CI-1040"],
       "ravoxertinib": ["GDC-0994"], "abemaciclib": ["LY2835219"], "ribociclib": ["LEE011"]}
CELL = {"NTERA2": "NTERA2CLD1", "NCIH3255": "H3255", "P31/FUJ": "P31FUJ"}
rows = []
for d in F.drug.unique():
    keys = {norm(d)} | {norm(s) for s in SYN.get(d, [])}
    m = L[L.drug_norm.isin(keys)]
    rows.append(dict(drug=d, gdsc1=";".join(sorted(set(m[m.source=="GDSC1"].drug_id))), gdsc2=";".join(sorted(set(m[m.source=="GDSC2"].drug_id))),
                     prism=";".join(sorted(set(m[m.source=="PRISM"].drug)))))
DM = pd.DataFrame(rows); DM.to_csv("/workspace/results/b_I23/tables/T0_drug_map.tsv", sep="\t", index=False)
out = []
for _, r in F[["drug", "cell"]].drop_duplicates().iterrows():
    keys = {norm(r.drug)} | {norm(s) for s in SYN.get(r.drug, [])}
    c = CELL.get(r.cell, r.cell)
    m = L[L.drug_norm.isin(keys) & (L.cell_norm == c)]
    for src, g in m.groupby("source"):
        out.append(dict(drug=r.drug, cell=r.cell, source=src, auc=g.auc.mean(), ln_ic50=g.ln_ic50.mean(), n_ids=len(g),
                        max_conc_uM=g.max_conc_uM.mean() if src != "PRISM" else np.nan))
S = pd.DataFrame(out); S.to_csv("/workspace/results/b_I23/tables/T2_sens_matched.tsv", sep="\t", index=False)
print(S.groupby("source").size()); print(DM.astype(bool).sum())
