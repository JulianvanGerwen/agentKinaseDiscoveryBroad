"""audit_D 05: LAMTOR1 S98 / RRAGC S15 (+ EPN1 S435 class view). Unit = study x agent x cell (dose>=100 nM; dose series only 1 uM).
Class medians, LOSO, genome-wide percentile of class median (sites with same n units), same-protein sites, S6Ki efficacy check."""
import common, numpy as np, pandas as pd
pd.set_option("display.width", 260); pd.set_option("display.max_rows", 400); pd.set_option("display.max_columns", 40)
OUT = "/workspace/results/audit_D/tables"
C = common.conds(); Z = common.zmat()
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
C["drug_class"] = np.where(C.single, C.agents.map(common.A2C), None)
D = C[C.usable & C.drug_class.notna() & ~(C.dose_nM < 100) & (C.dose_series_id.isna() | (C.dose_nM == 1000))]
D = D[~D.contrast_reversed.fillna(False).astype(bool) & D.background_agents.isna()]
print(D.groupby("drug_class").agg(n_cond=("conditionID", "size"), n_study=("study_pmid", "nunique")).to_string())
# unit matrix: sites x units (mean z over conditions in unit)
D["unit"] = D.drug_class + "|" + D.study_pmid + "|" + D.agents + "|" + D.cell_line_norm.fillna("NA")
U = Z[D.conditionID].T.groupby(D.set_index("conditionID").unit).mean().T
U.columns.name = None
uinfo = D.drop_duplicates("unit").set_index("unit")[["drug_class", "study_pmid", "agents", "cell_line_norm"]]
U.to_parquet("/workspace/cache/audit_D/unit_z.parquet"); uinfo.to_csv("/workspace/cache/audit_D/unit_info.tsv", sep="\t")
lam = [s for s in S.index if s.startswith("Q6IAA8_")]; rag = [s for s in S.index if s.startswith("Q9HB90_")]
epn = [s for s in S.index if s.startswith("Q9Y6I3_")]
PAN = {"Q6IAA8_S98": "LAMTOR1 S98", "Q9HB90_S15": "RRAGC S15", "Q9Y6I3_S435": "EPN1 S435", "Q96B36_S183": "PRAS40 S183",
       "Q13541_T46": "4EBP1 T46", "Q13541_S65": "4EBP1 S65", "Q6PKG0_S766": "LARP1 S766", "Q9HB90_S21": "RRAGC S21",
       "P62753_S240": "RPS6 S240", "P62753_S244": "RPS6 S244", "P23588_S422": "EIF4B S422", "Q96B36_T246": "PRAS40 T246",
       "P49841_S9": "GSK3B S9", "P27361_Y204": "ERK1 Y204", "P28482_Y187": "ERK2 Y187", "Q15418_S363": "RSK1 S363"}
for s in lam + rag + epn:
    PAN.setdefault(s, S.loc[s, "gene"] + " " + s.split("_")[1])
classes = ["cat_mTORi", "dual_PI3K_mTORi", "rapalog", "S6Ki", "AKTi", "PI3Ki", "MEKi", "ERKi", "RAFi", "EGFRi", "CDKi", "CK2i"]
rows = []
for cl in classes:
    cu = uinfo.index[uinfo.drug_class == cl]; Uc = U[cu]
    med = Uc.median(axis=1); n = Uc.notna().sum(axis=1)
    for s, lab in PAN.items():
        if s not in U.index or n[s] == 0: continue
        v = Uc.loc[s].dropna(); st = uinfo.loc[v.index, "study_pmid"]
        loso = [v[st != x].median() for x in st.unique()] if st.nunique() > 1 else [np.nan]
        # genome-wide percentile among sites measured in >= same number of units
        peer = med[n >= max(2, n[s])] if n[s] >= 2 else med[n >= 1]
        rows.append(dict(site=lab, cls=cl, n_units=int(n[s]), n_studies=st.nunique(), med_z=v.median(), n_lt_m2=int((v < -2).sum()),
                         loso_max=np.nanmax(loso), loso_min=np.nanmin(loso), pct_vs_sites=(peer < v.median()).mean() * 100, n_peer=len(peer)))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/class_effects_loso.tsv", sep="\t", index=False)
for cl in ["cat_mTORi", "dual_PI3K_mTORi", "rapalog", "S6Ki", "AKTi", "PI3Ki", "MEKi", "ERKi", "RAFi", "CDKi"]:
    print("\n==", cl); print(R[R.cls == cl].drop(columns="cls").round(2).to_string(index=False))
# per-unit table for cat_mTORi/dual/rapalog/S6Ki/AKTi/PI3Ki for candidates and controls
keys = ["Q6IAA8_S98", "Q6IAA8_S97", "Q9HB90_S15", "Q96B36_S183", "Q13541_T46", "P62753_S240", "P23588_S422", "Q96B36_T246"]
T = U.loc[[k for k in keys if k in U.index]].T.join(uinfo)
T = T[T.drug_class.isin(["cat_mTORi", "dual_PI3K_mTORi", "rapalog", "S6Ki", "AKTi", "PI3Ki"])]
T = T.rename(columns=PAN).dropna(subset=["LAMTOR1 S98", "RRAGC S15"], how="all").sort_values(["drug_class", "study_pmid"])
T.to_csv(f"{OUT}/units_mtor_candidates.tsv", sep="\t"); print(T.round(2).to_string())
