"""val_I03 step 4: CPTAC tumour co-variation (protein-adjusted, tumour-centred residuals; tumours only).
For each candidate: per-cohort Pearson r with pathway readouts and partial r (multiple regression on all
readouts of the relevant family); meta = Stouffer of Fisher z*sqrt(n-3).
Readouts: 4EBP1 T37/T46 (direct mTORC1), RPS6 S235/S236 (S6K), PRAS40 T246 (AKT), ERK1/2 pT/pY (MEK->ERK),
CK2 composite = mean of PSP CSNK2A1 substrate sites (>=5 measured), CDK composite (PSP CDK1 sites; mitosis proxy).
Out: tables/cptac_covariation.tsv, tables/cptac_meta.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/val_I03/tables"
CANDS = {"RRAGC S15": "Q9HB90_S15", "RRAGC S21": "Q9HB90_S21", "LAMTOR1 S98": "Q6IAA8_S98", "LAMTOR1 S97": "Q6IAA8_S97",
         "TBC1D5 S730": "Q92609_S730", "PRKAA1 S356": "Q13131_S356", "KPNA4 S60": "O00629_S60", "EPN1 S435": "Q9Y6I3_S435",
         # controls run through the same pipeline
         "4EBP1 S65": "Q13541_S65", "PRAS40 S183": "Q96B36_S183", "LARP1 S766": "Q6PKG0_S766", "ULK1 S758": "O75385_S758",
         "RPS6 S240": "P62753_S240", "EIF4B S422": "P23588_S422", "GSK3B S9": "P49841_S9", "TSC2 S939": "P49815_S939",
         "RSK1 S380": "Q15418_S380"}
READ = {"4EBP1_T37T46": ["Q13541_T37", "Q13541_T46"], "RPS6_S235S236": ["P62753_S235", "P62753_S236"],
        "PRAS40_T246": ["Q96B36_T246"], "ERK_pTpY": ["P27361_T202", "P27361_Y204", "P28482_T185", "P28482_Y187"]}
S = load.sites().set_index("site"); pk = S.psp_kinases_human.fillna("")
CK2set = set(S.index[pk.str.contains("CSNK2A1")]) - {"O00629_S60"}
CDKset = set(S.index[pk.str.contains("CDK1")]) - {"Q9Y6I3_S435"}
FAM = {"mTOR": ["4EBP1_T37T46", "RPS6_S235S236", "PRAS40_T246"], "MAPK": ["ERK_pTpY", "CK2_comp", "CDK_comp"]}

rows = []
for coh in cl.COHORTS:
    X = cl.matrix(coh, "input_tc_padj")
    X = X[[c for c in X.columns if c.endswith("_T")]]
    Xs = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
    R = {}
    for k, ss in READ.items():
        ss = [s for s in ss if s in Xs.index]
        if ss: R[k] = Xs.loc[ss].mean(0)
    for nm, st in [("CK2_comp", CK2set), ("CDK_comp", CDKset)]:
        ss = [s for s in st if s in Xs.index and Xs.loc[s].notna().mean() > .5]
        if len(ss) >= 5: R[nm] = Xs.loc[ss].mean(0)
    R = pd.DataFrame(R)
    for lab, s in CANDS.items():
        if s not in Xs.index: continue
        y = Xs.loc[s]
        for k in R.columns:
            d = pd.concat([y, R[k]], axis=1).dropna()
            if len(d) < 20: continue
            rows.append(dict(cohort=coh, site=lab, readout=k, kind="r", r=d.corr().iloc[0, 1], n=len(d)))
        for fam, ks in FAM.items():
            ks = [k for k in ks if k in R.columns]
            d = pd.concat([y.rename("y"), R[ks]], axis=1).dropna()
            if len(d) < 25 or len(ks) < 2: continue
            for k in ks:
                oth = [o for o in ks if o != k]
                A = np.c_[np.ones(len(d)), d[oth].values]
                ry = d.y.values - A @ np.linalg.lstsq(A, d.y.values, rcond=None)[0]
                rx = d[k].values - A @ np.linalg.lstsq(A, d[k].values, rcond=None)[0]
                rows.append(dict(cohort=coh, site=lab, readout=k, kind=f"partial_{fam}", r=np.corrcoef(ry, rx)[0, 1], n=len(d)))
    print(coh, R.shape)
T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/cptac_covariation.tsv", sep="\t", index=False)
T["zf"] = np.arctanh(T.r.clip(-.999, .999)) * np.sqrt(T.n - 3)
M = T.groupby(["site", "readout", "kind"]).agg(meta_z=("zf", "sum"), n_coh=("zf", "size"), r_med=("r", "median"),
                                                 n_pos=("r", lambda x: int((x > .1).sum()))).reset_index()
M["meta_z"] = M.meta_z / np.sqrt(M.n_coh)
M.to_csv(f"{OUT}/cptac_meta.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(M.pivot_table(index="site", columns=["kind", "readout"], values="meta_z").round(1).to_string())
print(M.pivot_table(index="site", columns=["kind", "readout"], values="r_med").round(2).to_string())
print(M[M.kind == "r"].pivot_table(index="site", columns="readout", values="n_coh").to_string())
