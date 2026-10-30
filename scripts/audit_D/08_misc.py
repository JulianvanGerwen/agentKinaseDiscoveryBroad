"""audit_D 08: misc checks - which units drive single-unit LAMTOR1 sites; RRAGC S15 in MEKi units; LAMTOR1 site windows/motifs;
AKTi/mTORi and PI3Ki/mTORi class-median ratios for candidates vs controls."""
import common, numpy as np, pandas as pd
pd.set_option("display.width", 250)
U = pd.read_parquet("/workspace/cache/audit_D/unit_z.parquet")
ui = pd.read_csv("/workspace/cache/audit_D/unit_info.tsv", sep="\t", index_col=0, dtype={"study_pmid": str})
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
lam = [s for s in U.index if s.startswith("Q6IAA8_")]
cu = ui.index[ui.drug_class.isin(["cat_mTORi", "dual_PI3K_mTORi"])]
print(U.loc[lam, cu].dropna(how="all", axis=1).dropna(how="all").round(2).T.to_string())
for s in ["Q6IAA8_S96", "Q6IAA8_S97", "Q6IAA8_S98", "Q6IAA8_T100"]:
    if s in S.index: print(s, S.loc[s, "sequence_window"], "n_cond", S.loc[s, "n_conditions"], "n_studies", S.loc[s, "n_studies"], "MTOR pct", M.loc[s, "MTOR"] if s in M.index else None)
mk = ui.index[ui.drug_class.isin(["MEKi"])]
print("RRAGC S15 MEKi units:", U.loc["Q9HB90_S15", mk].dropna().round(2).to_dict())
R = pd.read_csv("/workspace/results/audit_D/tables/class_effects_loso.tsv", sep="\t").pivot_table(index="site", columns="cls", values="med_z")
R["AKTi/cat"] = R.AKTi / R.cat_mTORi; R["PI3Ki/cat"] = R.PI3Ki / R.cat_mTORi; R["S6Ki/cat"] = R.S6Ki / R.cat_mTORi
print(R.loc[["LAMTOR1 S98", "RRAGC S15", "PRAS40 S183", "4EBP1 S65", "4EBP1 T46", "LARP1 S766", "EIF4B S422", "GSK3B S9", "PRAS40 T246"],
            ["cat_mTORi", "AKTi", "PI3Ki", "S6Ki", "AKTi/cat", "PI3Ki/cat", "S6Ki/cat"]].round(2).to_string())
