"""val_I03 step 3: decryptM dose-response (b_I17 curves_all) for candidates vs controls.
Key discriminator: pictilisib (PI3K-selective; inhibits mTOR only at ~uM) separates PI3K->AKT-driven sites
(high pEC50) from sites that need direct mTOR inhibition (low pEC50). AZD8055/dactolisib give potency ranks.
Out: tables/dose_pEC50.tsv (panel sites), tables/dose_reference_sets.tsv (PSP MTOR/AKT/S6K/ERK set medians)
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
src = open("/workspace/scripts/val_I03/01_panel.py").read()
exec(src[src.index("PANEL = {"):src.index("CLS = {")])
OUT = "/workspace/results/val_I03/tables"
d = pd.read_parquet("/workspace/results/b_I17/curves_all.parquet")
drugs = ["AZD8055", "dactolisib", "pictilisib", "MK2206", "PD0325901", "refametinib", "selumetinib",
         "afatinib", "gefitinib", "lapatinib", "dasatinib"]
d = d[d.agent.isin(drugs)]
P = d[d.site.isin(PANEL)].copy()
P["label"] = P.site.map(lambda s: PANEL[s][0]); P["group"] = P.site.map(lambda s: PANEL[s][1])
P["direction"] = np.sign(P.delta)
P = P[["exp_id", "agent", "cell_line", "time_min", "label", "group", "pEC50", "delta", "R2", "q", "sig", "ec_in_range", "max_abs"]]
P.to_csv(f"{OUT}/dose_pEC50.tsv", sep="\t", index=False)
# reference sets
import sys; sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
S = load.sites().set_index("site")
pk = S.psp_kinases_human.fillna("")
sets = {"PSP_MTOR": "MTOR", "PSP_AKT": "AKT1|AKT2|AKT3", "PSP_S6K": "RPS6KB1", "PSP_ERK": "MAPK1|MAPK3", "PSP_CK2": "CSNK2A1"}
rows = []
for nm, pat in sets.items():
    mem = set(S.index[pk.str.contains(pat)])
    q = d[d.site.isin(mem) & d.sig & d.ec_in_range & (d.delta < 0)]
    for (e, a), g in q.groupby(["exp_id", "agent"]):
        rows.append(dict(set=nm, exp_id=e, agent=a, n_sig_down=len(g), pEC50_med=g.pEC50.median(),
                         pEC50_q25=g.pEC50.quantile(.25), pEC50_q75=g.pEC50.quantile(.75)))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/dose_reference_sets.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
show = P[P.label.isin(["LAMTOR1 S97", "LAMTOR1 S98", "LAMTOR1 T100", "RRAGC S15", "TBC1D5 S730", "PRKAA1 S356",
                       "PRAS40 S183", "LARP1 S766", "RRAGC T394", "4EBP1 S65", "4EBP1 T46", "PRAS40 T246", "TSC2 S939",
                       "GSK3B S9", "RPS6 S244", "RPS6 S240", "EIF4B S422", "EPN1 S435", "KPNA4 S60", "ERK2 Y187",
                       "ERK1 Y204", "RSK1 S380", "RSK1 S363"]) & P.agent.isin(drugs)]
print(show.pivot_table(index=["group", "label"], columns=["agent", "exp_id"], values="pEC50").round(2).to_string())
print(show[show.sig].pivot_table(index=["group", "label"], columns=["agent", "exp_id"], values="pEC50").round(2).to_string())
print(R[R.agent.isin(["AZD8055", "dactolisib", "pictilisib", "MK2206"])].round(2).to_string())
