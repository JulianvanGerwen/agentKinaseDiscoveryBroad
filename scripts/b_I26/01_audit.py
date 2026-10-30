"""I26 step 1: main audit (all LoF studies). Edge classes + untestable reasons + random-relabel control.
Outputs: cache/b_I26/L_main.parquet (group,site,study,z,ncond), W_main.tsv (worked flags),
results tables/edge_audit.tsv, tables/worked_conditions.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *

E, Z, GC = load_base()
gold = gold_sets(E)
L, W = study_long(Z, GC, gold)
L["group"] = L.group.astype("category"); L["study"] = L.study.astype("category")
L.to_parquet(f"{C26}/L_main.parquet")
W.to_csv(f"{R26}/tables/worked_conditions.tsv", sep="\t", index=False, float_format="%.3g")
print(W.groupby("group").agg(n_cond=("worked", "size"), n_worked=("worked", "sum"),
                             st_worked=("study", lambda s: s[W.loc[s.index, "worked"]].nunique()),
                             n_st=("study", "nunique")).to_string())
A = classify(L)
# edge table: PSP edges of panel kinases
EP = E[E.group.notna()].copy()
EP = EP.merge(A, on=["group", "site"], how="left")
# untestable reasons
st_all = GC.groupby("group").study.nunique()
meas_any = {}
for g, d in GC.groupby("group"):
    zz = Z[d.conditionID.unique()]
    stc = d.drop_duplicates("conditionID").set_index("conditionID").study
    m = zz.notna().T.groupby(stc.reindex(zz.columns).values).any().sum()   # n studies measuring each site
    meas_any[g] = m
EP["n_lof_studies"] = EP.group.map(st_all)
EP["n_st_measured_any"] = [meas_any[g].get(s, 0) for g, s in zip(EP.group, EP.site)]
EP["cls"] = EP.cls.fillna("untestable")
EP["n_meas"] = EP.n_meas.fillna(0).astype(int)
EP["reason"] = np.select([EP.n_lof_studies < 2, EP.n_st_measured_any < 2, EP.cls.eq("untestable")],
                         ["kinase<2 LoF studies", "site measured <2 LoF studies", "measured but inhibition not verified in >=2 studies"], "")
EP["is_gold"] = [s in gold.get(g, ()) for g, s in zip(EP.group, EP.site)]
EP.to_csv(f"{R26}/tables/edge_audit.tsv", sep="\t", index=False, float_format="%.3g")
print(EP.cls.value_counts().to_dict())
print(EP.reason.value_counts().to_dict())
print(EP[EP.cls != "untestable"].groupby("group").cls.value_counts().unstack(fill_value=0).to_string())
