"""val_I03 step 5: T2 (CK2 -> KPNA4 S60) and T3 (ERK -> EPN1 S435) specificity tests.
 - percentile of the candidate among all measured sites, and among PSP CK2 / ERK substrates, in CK2i / MEKi / ERKi conditions
 - 'pathway-verified' MEKi/ERKi/RAFi/EGFRi units: does the candidate fall when ERK pT/pY (MEKi/RAFi/EGFRi) or
   RSK1 S363 (ERKi) clearly falls? (unit-level concordance)
 - motif check (CK2 acidic: D/E at +1..+3; ERK: S-P with P at -2)
 - EGF/FGF stimulation effects
Out: tables/t2_t3_percentiles.tsv, tables/t3_mapk_concordance.tsv, tables/motif_check.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
CA = "/workspace/cache/val_I03"; OUT = "/workspace/results/val_I03/tables"
C = pd.read_csv(f"{CA}/cond_class.tsv", sep="\t", low_memory=False, dtype={"study_pmid": str})
C = C[C.usable & C.drug_class.isin(["CK2i", "MEKi", "ERKi", "RAFi", "EGFRi", "cat_mTORi", "dual_PI3K_mTORi", "AKTi", "S6Ki"]) &
      ~(C.dose_nM < 100) & (C.dose_series_id.isna() | (C.dose_nM == 1000))]
X = load.sites_matrix(conditions=list(C.cid))
med = X.median(); mad = (X - med).abs().median() * 1.4826; Z = (X - med) / mad
S = load.sites().set_index("site"); pk = S.psp_kinases_human.fillna("")
ck2 = S.index[pk.str.contains("CSNK2A1")].intersection(Z.index)
erk = S.index[pk.str.contains("MAPK1|MAPK3")].intersection(Z.index)
mtor = S.index[pk.str.contains("MTOR")].intersection(Z.index)
CAND = {"KPNA4 S60": "O00629_S60", "EPN1 S435": "Q9Y6I3_S435", "LAMTOR1 S98": "Q6IAA8_S98", "LAMTOR1 S97": "Q6IAA8_S97",
        "RRAGC S15": "Q9HB90_S15", "TBC1D5 S730": "Q92609_S730", "RSK1 S363": "Q15418_S363"}
rows = []
for cid in C.cid:
    z = Z[cid].dropna()
    if len(z) < 500: continue
    for lab, s in CAND.items():
        if s not in z.index: continue
        v = z[s]
        rows.append(dict(cid=cid, site=lab, z=v, pct_all=(z < v).mean() * 100,
                         pct_CK2=(z.reindex(ck2).dropna() < v).mean() * 100,
                         pct_ERK=(z.reindex(erk).dropna() < v).mean() * 100,
                         pct_MTOR=(z.reindex(mtor).dropna() < v).mean() * 100,
                         CK2set_med=z.reindex(ck2).median(), ERKset_med=z.reindex(erk).median()))
P = pd.DataFrame(rows).merge(C[["cid", "study_pmid", "agents", "drug_class", "cell_line_norm", "time_min"]], on="cid")
P.to_csv(f"{OUT}/t2_t3_percentiles.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(P[(P.site == "KPNA4 S60") & (P.drug_class == "CK2i")].round(2).to_string())
print(P[P.site.isin(["KPNA4 S60", "EPN1 S435"])].groupby(["site", "drug_class"])[["z", "pct_all"]].median().round(2))

# T3 concordance: per unit, ERK loop drop vs EPN1 S435 drop
L = pd.read_parquet(f"{CA}/panel_long.parquet")
U = pd.read_csv(f"{OUT}/unit_effects.tsv", sep="\t", dtype={"study_pmid": str})
W = U.pivot_table(index=["drug_class", "study_pmid", "agents", "cell_line_norm"], columns="label", values="z")
loop = W[[c for c in ["ERK1 T202", "ERK1 Y204", "ERK2 T185", "ERK2 Y187"] if c in W]].mean(1)
T = pd.DataFrame({"ERK_loop": loop, "RSK1_S363": W.get("RSK1 S363"), "RSK1_S380": W.get("RSK1 S380"),
                  "EPN1_S435": W["EPN1 S435"], "KPNA4_S60": W["KPNA4 S60"], "EPN1_S454": W.get("EPN1 S454")}).reset_index()
T = T[T.drug_class.isin(["MEKi", "ERKi", "RAFi", "EGFRi"])]
T["pathway_down"] = np.where(T.drug_class == "ERKi", T.RSK1_S363 < -2, T.ERK_loop < -2)
T.to_csv(f"{OUT}/t3_mapk_concordance.tsv", sep="\t", index=False)
print(T.round(2).to_string())
for nm in ["EPN1_S435", "KPNA4_S60", "EPN1_S454"]:
    a = T[T.pathway_down][nm].dropna(); b = T[~T.pathway_down][nm].dropna()
    print(nm, "pathway-down units: median z %.2f (n=%d, %d<-2) | not-down: %.2f (n=%d, %d<-2)" %
          (a.median(), len(a), (a < -2).sum(), b.median(), len(b), (b < -2).sum()))
    d = T[["ERK_loop", nm]].dropna(); print("   spearman ERK_loop vs", nm, round(d.corr("spearman").iloc[0, 1], 2), len(d))

# motif check
def motif(w):
    c = w[7]; r = {"window": w, "res": c}
    r["CK2_acidic_+1..+3"] = sum(x in "DE" for x in w[8:11]); r["CK2_+3_acidic"] = w[10] in "DE"
    r["Pro+1"] = w[8] == "P"; r["Pro-2"] = w[5] == "P"; r["Arg-3"] = w[4] == "R"; r["Arg-2"] = w[5] == "R"
    r["hydrophobic+1"] = w[8] in "FLIMVWYP"
    return r
M = pd.DataFrame({lab: motif(S.loc[s, "sequence_window"]) for lab, s in CAND.items()}).T
M.to_csv(f"{OUT}/motif_check.tsv", sep="\t"); print(M)
