"""audit_D 06: decryptM curves for LAMTOR1/RRAGC/EPN1/KPNA4 vs controls; S6K-inhibitor efficacy check; motif pct."""
import common, numpy as np, pandas as pd
pd.set_option("display.width", 260); pd.set_option("display.max_rows", 300)
OUT = "/workspace/results/audit_D/tables"
d = pd.read_parquet("/workspace/results/b_I17/curves_all.parquet")
P = {"Q6IAA8_S98": "LAMTOR1 S98", "Q6IAA8_S97": "LAMTOR1 S97", "Q6IAA8_S27": "LAMTOR1 S27", "Q9HB90_S15": "RRAGC S15",
     "Q96B36_S183": "PRAS40 S183", "Q96B36_T246": "PRAS40 T246", "Q13541_S65": "4EBP1 S65", "Q13541_T46": "4EBP1 T46",
     "Q6PKG0_S766": "LARP1 S766", "P23588_S422": "EIF4B S422", "P62753_S240": "RPS6 S240", "P62753_S244": "RPS6 S244", "P49841_S9": "GSK3B S9",
     "Q9Y6I3_S435": "EPN1 S435", "Q9Y6I3_S454": "EPN1 S454", "P27361_Y204": "ERK1 Y204", "P28482_Y187": "ERK2 Y187",
     "Q15418_S380": "RSK1 S380", "Q15418_S363": "RSK1 S363", "O00629_S60": "KPNA4 S60"}
q = d[d.site.isin(P) & d.agent.isin(["AZD8055", "dactolisib", "pictilisib", "MK2206", "PD0325901", "refametinib", "selumetinib", "lapatinib", "gefitinib", "afatinib", "dasatinib"])].copy()
q["label"] = q.site.map(P)
q["tag"] = q.agent + "|" + q.cell_line + "|" + q.time_min.astype(int).astype(str)
cols = ["pEC50", "delta", "R2", "sig"]
T = q.pivot_table(index="label", columns="tag", values="pEC50")
Dl = q.pivot_table(index="label", columns="tag", values="delta")
Sg = q.pivot_table(index="label", columns="tag", values="sig", aggfunc="max")
keep = [c for c in T.columns if c.split("|")[0] in ["AZD8055", "dactolisib", "pictilisib", "MK2206"]]
out = T[keep].round(2).astype(str) + " (" + Dl[keep].round(2).astype(str) + ("" ) + Sg[keep].map(lambda x: "*" if x is True or x == 1 else "").astype(str) + ")"
print("pEC50 (delta, *=sig) mTOR/PI3K/AKT drugs:"); print(out.to_string())
keep2 = [c for c in T.columns if c.split("|")[0] not in ["AZD8055", "dactolisib", "pictilisib", "MK2206"]]
lab2 = ["EPN1 S435", "EPN1 S454", "ERK1 Y204", "ERK2 Y187", "RSK1 S380", "RSK1 S363", "KPNA4 S60", "LAMTOR1 S98", "PRAS40 S183"]
out2 = T.loc[T.index.intersection(lab2), keep2].round(2).astype(str) + " (" + Dl.loc[T.index.intersection(lab2), keep2].round(2).astype(str) + Sg.loc[T.index.intersection(lab2), keep2].map(lambda x: "*" if x is True or x == 1 else "").astype(str) + ")"
print("\npEC50 MAPK-pathway drugs:"); print(out2.T.to_string())
q.to_csv(f"{OUT}/decryptM_panel.tsv", sep="\t", index=False)
# S6Ki efficacy: S6K substrate sites in S6Ki conditions
Z = common.zmat(); C = common.conds().set_index("conditionID")
s6 = C.index[C.agents.isin(["PF4708671", "LY2584702"])]
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
pk = S.psp_kinases_human.fillna("")
s6set = S.index[pk.str.contains("RPS6KB1")].intersection(Z.index)
rows = []
for c in s6:
    z = Z[c]; v = z.reindex(s6set).dropna()
    rows.append(dict(cid=c, agent=C.loc[c, "agents"], cell=C.loc[c, "cell_line_norm"], study=C.loc[c, "study_pmid"],
                     PSP_S6K_med=v.median(), n=len(v), RPS6_S235=z.get("P62753_S235"), RPS6_S236=z.get("P62753_S236"),
                     RPS6_S240=z.get("P62753_S240"), RPS6_S244=z.get("P62753_S244"), EIF4B_S422=z.get("P23588_S422"),
                     IRS1_S1101=z.get("P35568_S1101"), LAMTOR1_S98=z.get("Q6IAA8_S98")))
R = pd.DataFrame(rows); print("\nS6K inhibitor efficacy (z):"); print(R.round(2).to_string()); R.to_csv(f"{OUT}/s6ki_efficacy.tsv", sep="\t", index=False)
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
for s in ["Q6IAA8_S98", "Q9HB90_S15", "Q9Y6I3_S435"]:
    m = M.loc[s]; print("\n", P[s], S.loc[s, "sequence_window"], "top:", m.sort_values(ascending=False).head(12).to_dict())
    print("   ", {k: int(m[k]) for k in ["MTOR", "AKT1", "P70S6K", "RPS6KB1", "RPS6KA1", "RPS6KA3", "SGK1", "PRKACA", "ERK1", "MAPK1", "MAPK3", "CDK1", "CDK2", "GSK3B", "CSNK2A1"] if k in m.index})
