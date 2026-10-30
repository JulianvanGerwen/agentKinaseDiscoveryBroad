"""val_I03 step 1: extract candidate + control sites from the compendium, robust z per condition,
annotate conditions with drug classes. Also PSP kinase-set medians per condition.
Out: cache/val_I03/panel_long.parquet, set_long.parquet, cond_class.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/cache/val_I03"

PANEL = {
 # candidates
 "Q9HB90_S15": ("RRAGC S15", "cand_mTOR"), "Q6IAA8_S98": ("LAMTOR1 S98", "cand_mTOR"),
 "Q92609_S730": ("TBC1D5 S730", "cand_mTOR"), "Q13131_S356": ("PRKAA1 S356", "posctrl_mTOR"),
 "O00629_S60": ("KPNA4 S60", "cand_CK2"), "Q9Y6I3_S435": ("EPN1 S435", "cand_ERK"),
 # neighbours of candidates
 "Q6IAA8_S97": ("LAMTOR1 S97", "neighbour"), "Q6IAA8_T100": ("LAMTOR1 T100", "neighbour"),
 "Q92609_S732": ("TBC1D5 S732", "neighbour"), "Q9HB90_T10": ("RRAGC T10", "neighbour"),
 "Q9HB90_Y16": ("RRAGC Y16", "neighbour"), "Q9HB90_S95": ("RRAGC S95", "neighbour"),
 # direct mTORC1
 "Q13541_T37": ("4EBP1 T37", "direct_mTORC1"), "Q13541_T46": ("4EBP1 T46", "direct_mTORC1"),
 "Q13541_S65": ("4EBP1 S65", "direct_mTORC1"), "Q13541_T70": ("4EBP1 T70", "direct_mTORC1"),
 "O75385_S758": ("ULK1 S758", "direct_mTORC1"), "Q6PKG0_S766": ("LARP1 S766", "direct_mTORC1"),
 "Q6PKG0_S774": ("LARP1 S774", "direct_mTORC1"), "P19484_S122": ("TFEB S122", "direct_mTORC1"),
 "Q9HB90_S21": ("RRAGC S21", "direct_mTORC1"), "Q9HB90_T394": ("RRAGC T394", "direct_mTORC1"),
 "Q96B36_S183": ("PRAS40 S183", "direct_mTORC1"),
 # S6K
 "P62753_S235": ("RPS6 S235", "S6K"), "P62753_S236": ("RPS6 S236", "S6K"),
 "P62753_S240": ("RPS6 S240", "S6K"), "P62753_S244": ("RPS6 S244", "S6K"),
 "P23588_S422": ("EIF4B S422", "S6K"), "P35568_S1101": ("IRS1 S1101", "S6K"),
 # AKT
 "Q96B36_T246": ("PRAS40 T246", "AKT"), "P49815_T1462": ("TSC2 T1462", "AKT"),
 "P49815_S939": ("TSC2 S939", "AKT"), "P49815_S981": ("TSC2 S981", "AKT"),
 "P49841_S9": ("GSK3B S9", "AKT"), "P49840_S21": ("GSK3A S21", "AKT"),
 # mTORC2
 "P31749_S473": ("AKT1 S473", "mTORC2"), "Q92597_T346": ("NDRG1 T346", "mTORC2/SGK"),
 # ERK pathway
 "P28482_T185": ("ERK2 T185", "ERK_loop"), "P28482_Y187": ("ERK2 Y187", "ERK_loop"),
 "P27361_T202": ("ERK1 T202", "ERK_loop"), "P27361_Y204": ("ERK1 Y204", "ERK_loop"),
 "Q15418_S380": ("RSK1 S380", "ERK_down"), "Q15418_S363": ("RSK1 S363", "ERK_direct"),
 # EPN1 other sites
 "Q9Y6I3_S454": ("EPN1 S454", "neighbour"), "Q9Y6I3_T470": ("EPN1 T470", "neighbour"),
 "Q9Y6I3_T460": ("EPN1 T460", "neighbour"),
}

CLS = {
 "cat_mTORi": "AZD8055 KU0063794 torin1 torin2 INK128 sapanisertib PP242 torkinib mTORi",
 "dual_PI3K_mTORi": "dactolisib omipalisib PI103",
 "rapalog": "rapamycin everolimus temsirolimus",
 "S6Ki": "PF4708671 LY2584702",
 "AKTi": "MK2206 capivasertib ipatasertib GSK690693 akti_akt12 aktVIII AT13148 perifosine",
 "PI3Ki": "pictilisib AZD6482 idelalisib TGX221 HS173 PIK294 PI3Ka_inh2 wortmannin LY294002 PI3Ki",
 "CK2i": "CX4945 TBCA",
 "MEKi": "selumetinib PD0325901 trametinib cobimetinib refametinib U0126 MEKi PD184352",
 "ERKi": "ERKi_SCH772984 ERKi_generic ravoxertinib ulixertinib",
 "RAFi": "vemurafenib dabrafenib SB590885 BRAFi LY3009120",
 "EGFRi": "gefitinib erlotinib lapatinib afatinib AG1478 AZD3759 PD153035 PD168393 CUDC101",
 "CDKi": "RO3306 dinaciclib AZD5438 flavopiridol",
 "p38i": "PH797804 SB202190 SB203580 SB239063 TAK715 p38i BIRB796",
 "GSK3i": "CHIR99021 GSK3i_BIO LY2090314 SB216763 tideglusib",
}
A2C = {a: c for c, s in CLS.items() for a in s.split()}

C = load.conditions()
C["cid"] = C.conditionID
C["drug_class"] = pd.Series([None]*len(C), dtype=object)
single = (~C.is_combination.astype(bool)) & C.agents.notna() & ~C.agents.fillna("").str.contains(";")
C.loc[single, "drug_class"] = C.loc[single, "agents"].map(A2C)
# stimulation classes
st = C.stimulus.fillna("")
C.loc[C.drug_class.isna() & st.isin(["insulin"]), "drug_class"] = "stim_insulin"
C.loc[C.drug_class.isna() & st.isin(["serum"]), "drug_class"] = "stim_serum"
C.loc[C.drug_class.isna() & st.str.contains("starvation"), "drug_class"] = "starvation"
C.loc[C.drug_class.isna() & st.str.contains("EGF|FGF|PMA|NGF"), "drug_class"] = "stim_GF_" + st
C["usable"] = ~C.qc_flag.isin(["duplicate_vector", "scale_suspect"])
C.to_csv(f"{OUT}/cond_class.tsv", sep="\t", index=False)
print(C.drug_class.value_counts().head(30))

X = load.sites_matrix()
print("matrix", X.shape)
V = X.values.astype(np.float32)
med = np.nanmedian(V, 0)
mad = np.nanmedian(np.abs(V - med), 0) * 1.4826
Z = (V - med) / mad
Zdf = pd.DataFrame(Z, index=X.index, columns=X.columns)

sel = [s for s in PANEL if s in X.index]
print("missing", [s for s in PANEL if s not in X.index])
L = X.loc[sel].stack().rename("log2fc").to_frame().join(Zdf.loc[sel].stack().rename("z"))
L.index.names = ["site", "cid"]; L = L.reset_index()
L["label"] = L.site.map(lambda s: PANEL[s][0]); L["group"] = L.site.map(lambda s: PANEL[s][1])
L.to_parquet(f"{OUT}/panel_long.parquet"); print(L.shape)

# PSP kinase sets (all sites with PSP-annotated kinase)
S = load.sites().set_index("site")
sets = {"PSP_MTOR": ["MTOR"], "PSP_S6K": ["RPS6KB1", "RPS6KB2"], "PSP_AKT": ["AKT1", "AKT2", "AKT3"],
        "PSP_CK2": ["CSNK2A1", "CSNK2A2"], "PSP_ERK": ["MAPK1", "MAPK3"], "PSP_CDK1": ["CDK1"],
        "PSP_RSK": ["RPS6KA1", "RPS6KA3"]}
pk = S.psp_kinases_human.fillna("").str.split(";")
rows = []
for nm, ks in sets.items():
    mem = [s for s in S.index[pk.map(lambda x: any(k in x for k in ks))] if s in Zdf.index and s not in PANEL]
    if nm == "PSP_ERK":  # proline-directed only
        mem = [s for s in mem if S.loc[s, "sequence_window"][8:9] == "P"]
    zz = Zdf.loc[mem]
    rows.append(pd.DataFrame({"set": nm, "cid": Zdf.columns, "z_med": zz.median(0).values,
                              "n": zz.notna().sum(0).values}))
    print(nm, len(mem))
pd.concat(rows).to_parquet(f"{OUT}/set_long.parquet")
