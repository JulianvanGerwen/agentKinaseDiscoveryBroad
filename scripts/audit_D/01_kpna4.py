"""audit_D 01: KPNA4 S60 vs CK2 inhibition. Known CK2 sites in the same CX-4945/TBCA conditions, rank of KPNA4 S60,
Cantley-motif-set behaviour (CK2 and Klaeger off-targets), other inhibitors hitting CX-4945 off-targets, dynamics."""
import common, numpy as np, pandas as pd, re
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500); pd.set_option("display.max_columns", 40)
OUT = "/workspace/results/audit_D/tables"
C = common.conds(); Z = common.zmat()
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
K = "O00629_S60"
ck2c = C[C.agents.isin(["CX4945", "TBCA"])].conditionID.tolist()
lf = common.lfc(conds=ck2c)
# --- (i) KPNA4 S60 rank + KPNA4 other sites
rows = []
for c in ck2c:
    z = Z[c].dropna(); a = C.set_index("conditionID").loc[c]
    r = dict(cid=c, agent=a.agents, cell=a.cell_line_norm, time=a.time_min, n_sites=len(z))
    if K in z.index:
        r.update(lfc=lf.loc[K, c], z=z[K], pct=(z < z[K]).mean() * 100, rank=int((z < z[K]).sum()) + 1)
    rows.append(r)
R = pd.DataFrame(rows); print(R.round(3).to_string()); R.to_csv(f"{OUT}/kpna4_ck2i_rank.tsv", sep="\t", index=False)
kp = [s for s in Z.index if s.startswith("O00629_")]
print("KPNA4 sites in CK2i conditions (z):"); print(Z.loc[kp, ck2c].dropna(how="all").round(2).to_string())
print(S.loc[[s for s in kp if s in S.index and Z.loc[s, ck2c].notna().any()], "sequence_window"])

# --- PSP CK2 substrates (human), LT evidence
ks = pd.read_csv("/workspace/data/Kinase_Substrate_Dataset", sep="\t", skiprows=3, low_memory=False, encoding="latin-1")
ks = ks[(ks.KIN_ORGANISM == "human") & (ks.SUB_ORGANISM == "human") & ks.GENE.isin(["CSNK2A1", "CSNK2A2"])]
ks["site"] = ks.SUB_ACC_ID + "_" + ks.SUB_MOD_RSD
psd = pd.read_csv("/workspace/data/Phosphorylation_site_dataset", sep="\t", skiprows=3, low_memory=False, encoding="latin-1")
psd = psd[psd.ORGANISM == "human"]; psd["site"] = psd.ACC_ID + "_" + psd.MOD_RSD.str.replace("-p", "")
lt = psd.groupby("site").LT_LIT.max()
ks["LT_LIT"] = ks.site.map(lt); ks["in_vivo"] = ks.IN_VIVO_RXN.astype(str).str.strip() == "X"
ck2 = ks.groupby("site").agg(gene=("SUB_GENE", "first"), LT=("LT_LIT", "max"), in_vivo=("in_vivo", "max")).reset_index()
ck2 = ck2[ck2.site.isin(Z.index)]
cx = [c for c in ck2c if C.set_index("conditionID").loc[c, "agents"] == "CX4945"]
T = ck2.set_index("site").join(Z.loc[ck2.site, ck2c])
T["n_meas_cx"] = T[cx].notna().sum(1)
Tm = T[T.n_meas_cx > 0].sort_values("LT", ascending=False)
Tm.to_csv(f"{OUT}/ck2_psp_sites_in_ck2i.tsv", sep="\t")
print("\nPSP CK2 substrates measured in >=1 CX-4945 condition:", len(Tm))
print(Tm.round(2).to_string())
summ = []
for c in ck2c:
    z = Z[c].dropna()
    for nm, mem in [("PSP_CK2_all", ck2.site), ("PSP_CK2_vivo_LT>=3", ck2.site[(ck2.LT >= 3) & ck2.in_vivo])]:
        v = z.reindex(mem).dropna()
        summ.append(dict(cid=c, set=nm, n=len(v), med_z=v.median(), frac_lt_m1=(v < -1).mean(), frac_lt_m2=(v < -2).mean(),
                         bg_frac_lt_m2=(z < -2).mean()))
SU = pd.DataFrame(summ); print(SU.round(3).to_string()); SU.to_csv(f"{OUT}/ck2_set_summary.tsv", sep="\t", index=False)

# --- (ii) Cantley motif sets
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
print("\nKPNA4 S60 top motif percentiles:", M.loc[K].sort_values(ascending=False).head(25).to_dict())
for k in ["CK2A1", "CSNK2A1", "CSNK2A2", "CK1A", "CSNK1A1", "CSNK1D", "CLK1", "CLK3", "DYRK1A", "MAPKAPK3", "PIM1", "GSK3B", "HIPK2", "PLK1", "CDK5", "SLK", "STK10", "MAPK1", "MAPK3", "GRK2"]:
    if k in M.columns: print(k, M.loc[K, k], end="; ")
print()
rows = []
for c in ck2c:
    z = Z[c].dropna(); zz = z.reindex(M.index); mm = M.loc[zz.dropna().index]; zz = zz.dropna()
    base = zz.median()
    for k in M.columns:
        top = zz[mm[k].values >= 98]
        rows.append(dict(cid=c, kinase=k, n=len(top), med_z=top.median() - base, frac_lt_m2=(top < -2).mean()))
MS = pd.DataFrame(rows)
MS.to_csv(f"{OUT}/motif_set_shift_ck2i.tsv", sep="\t", index=False)
W = MS.pivot_table(index="kinase", columns="cid", values="med_z")
W["mean_CX"] = W[cx].mean(1); W = W.sort_values("mean_CX")
print("\nMotif-top2% set median-z shift, most negative under CX-4945:"); print(W.head(25).round(2).to_string())
offt = ["CSNK2A1", "CSNK2A2", "CLK1", "CLK2", "CLK3", "CLK4", "DYRK1A", "DYRK1B", "DYRK2", "HIPK1", "HIPK2", "HIPK3", "PIM1", "PIM3",
        "GSK3A", "GSK3B", "MAPKAPK3", "STK10", "SLK", "CDK5", "PLK1", "CDK2", "CK1A", "CSNK1A1", "CSNK1D", "MAPK1"]
print(W.loc[[k for k in offt if k in W.index]].round(2).to_string())
print("rank of KPNA4-style: CK2 rank in list", list(W.index).index("CSNK2A1") if "CSNK2A1" in W.index else None, "of", len(W))
