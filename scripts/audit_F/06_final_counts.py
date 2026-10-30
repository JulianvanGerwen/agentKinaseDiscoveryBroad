"""Final per-site lab counts (verified inhibition, own markers), Thibault (HCT116 vemurafenib = paradoxical-activation context) excluded
from MAPK evidence; PKC and CDK separately."""
import pandas as pd, numpy as np
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
L = pd.read_csv("/workspace/results/audit_F/tables/site_condition_long.tsv", sep="\t", dtype={"study": str})
L = L[L.verified]
L = L[~((L.lab == "Thibault") & L.cls.isin(["RAF"]))]
grp = {"AKT": ["AKT"], "PKC": ["PKC"], "CDK": ["CDK"], "MTOR": ["MTOR"], "MEK/RAF": ["MEK", "RAF"], "ERKi": ["ERK"], "ALK": ["ALK"], "EGFR": ["EGFR"]}
rows = []
for (g, s), d in L.groupby(["gene", "site"]):
    r = dict(gene=g, site=s)
    for k, cl in grp.items():
        dd = d[d.cls.isin(cl)]
        st = dd.groupby(["lab", "study"]).z.mean().groupby("lab").mean()
        r[k] = f"{(st <= -1).sum()}dn/{(st > -0.5).sum()}flat/{len(st)}" if len(st) else "-"
        r[k + "_z"] = ";".join(f"{l[:4]}{v:+.1f}" for l, v in st.items())
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv("/workspace/results/audit_F/tables/site_lab_counts_final.tsv", sep="\t", index=False)
print(R[["gene"] + list(grp)].to_string())
print(R[["gene"] + [k + "_z" for k in ["AKT", "PKC", "CDK", "MEK/RAF", "ERKi", "ALK", "MTOR"]]].to_string())
