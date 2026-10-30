"""val_I03 step 2: compendium discriminating tests (T1 a/b/c/e, T2, T3).
Inputs: cache/val_I03/panel_long.parquet (robust per-condition z), cond_class.tsv
Outputs (results/val_I03/tables):
  class_effects.tsv        site x drug class: median z over study-agent-cell units, n units, n units z<-2
  matched_ratios.tsv       within-study matched contrasts (rapalog/catalytic, S6Ki/mTORi, AKTi/mTORi, PI3Ki/mTORi...)
  profile_corr.tsv         compendium-wide profile correlations of candidates with composite readouts (+partial)
  mk2206_timing.tsv        A549 MK2206 30 vs 60 min (same doses) for candidates vs controls
  stim_effects.tsv         stimulation conditions (insulin/serum/GF/starvation)
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
CA = "/workspace/cache/val_I03"; OUT = "/workspace/results/val_I03/tables"
L = pd.read_parquet(f"{CA}/panel_long.parquet")
C = pd.read_csv(f"{CA}/cond_class.tsv", sep="\t", low_memory=False, dtype={"study_pmid": str})
L = L.merge(C[["cid", "study_pmid", "agents", "drug_class", "cell_line_norm", "dose_nM", "time_min",
               "dose_series_id", "time_series_id", "stimulus", "usable", "qc_flag"]], on="cid")
L = L[L.usable]
CANDS = ["RRAGC S15", "LAMTOR1 S98", "TBC1D5 S730", "PRKAA1 S356", "KPNA4 S60", "EPN1 S435"]

# ---------- class effects (units = study x agent x cell; drop low doses <100 nM and dose-series except 1 uM) ----------
D = L[L.drug_class.notna() & ~L.drug_class.str.startswith("stim") & (L.drug_class != "starvation")].copy()
D = D[~(D.dose_nM < 100)]
D = D[D.dose_series_id.isna() | (D.dose_nM == 1000)]
U = D.groupby(["label", "group", "drug_class", "study_pmid", "agents", "cell_line_norm"]).z.mean().reset_index()
U.to_csv(f"{OUT}/unit_effects.tsv", sep="\t", index=False)
CE = U.groupby(["label", "group", "drug_class"]).z.agg(z_med="median", n_units="size",
                                                         n_down=lambda x: int((x < -2).sum())).reset_index()
CE.to_csv(f"{OUT}/class_effects.tsv", sep="\t", index=False)
wide = CE.pivot_table(index=["group", "label"], columns="drug_class", values="z_med").round(2)
print(wide[["cat_mTORi", "dual_PI3K_mTORi", "rapalog", "S6Ki", "AKTi", "PI3Ki", "MEKi", "ERKi", "RAFi", "CK2i", "CDKi"]].to_string())

# ---------- matched within-study contrasts ----------
def matched(num, den, sub=None):
    rows = []
    a = U[U.drug_class == num].groupby(["label", "group", "study_pmid", "cell_line_norm"]).z.mean()
    b = U[U.drug_class == den].groupby(["label", "group", "study_pmid", "cell_line_norm"]).z.mean()
    j = pd.concat([a.rename("z_num"), b.rename("z_den")], axis=1).dropna().reset_index()
    j["contrast"] = f"{num}/{den}"
    return j
M = pd.concat([matched("rapalog", "cat_mTORi"), matched("S6Ki", "cat_mTORi"), matched("AKTi", "cat_mTORi"),
               matched("PI3Ki", "cat_mTORi"), matched("dual_PI3K_mTORi", "cat_mTORi"),
               matched("MEKi", "cat_mTORi"), matched("CK2i", "MEKi"), matched("CK2i", "ERKi"),
               matched("ERKi", "MEKi"), matched("S6Ki", "ERKi")])
# ratio is only interpretable when the denominator effect is clear
M["ratio"] = np.where(M.z_den < -2, M.z_num / M.z_den, np.nan)
M.to_csv(f"{OUT}/matched_ratios.tsv", sep="\t", index=False)
print(M[M.label.isin(CANDS + ["4EBP1 T37", "4EBP1 T46", "4EBP1 S65", "PRAS40 S183", "ULK1 S758", "LARP1 S766",
                              "RPS6 S235", "RPS6 S236", "RPS6 S240", "RPS6 S244", "EIF4B S422", "PRAS40 T246",
                              "GSK3B S9"]) & M.contrast.isin(["rapalog/cat_mTORi", "S6Ki/cat_mTORi", "AKTi/cat_mTORi"])]
      .sort_values(["contrast", "group", "label"]).round(2).to_string())

# ---------- compendium-wide profile correlations ----------
W = L.pivot_table(index="cid", columns="label", values="z")
W = W.clip(-10, 10)
def comp(labels):
    return W[[l for l in labels if l in W]].mean(1)
R = {"direct_4EBP1": comp(["4EBP1 T37", "4EBP1 T46"]),
     "direct_mTORC1": comp(["4EBP1 T37", "4EBP1 T46", "4EBP1 S65", "PRAS40 S183", "LARP1 S766", "ULK1 S758"]),
     "S6K": comp(["RPS6 S235", "RPS6 S236", "RPS6 S240", "RPS6 S244", "EIF4B S422"]),
     "AKT": comp(["PRAS40 T246", "TSC2 T1462", "TSC2 S939", "GSK3B S9", "GSK3A S21"]),
     "ERK_loop": comp(["ERK1 T202", "ERK1 Y204", "ERK2 T185", "ERK2 Y187"]),
     "RSK_S380": comp(["RSK1 S380"])}
S = pd.read_parquet(f"{CA}/set_long.parquet").pivot_table(index="cid", columns="set", values="z_med")
R["PSP_CK2_set"] = S["PSP_CK2"].reindex(W.index); R["PSP_ERK_set"] = S["PSP_ERK"].reindex(W.index)
R["PSP_CDK1_set"] = S["PSP_CDK1"].reindex(W.index)
R = pd.DataFrame(R)

def pcor(y, x, zs):
    d = pd.concat([y, x] + zs, axis=1).dropna()
    if len(d) < 15: return np.nan, len(d)
    r = np.linalg.lstsq(np.c_[np.ones(len(d)), d.iloc[:, 2:].values], d.iloc[:, :2].values, rcond=None)
    res = d.iloc[:, :2].values - np.c_[np.ones(len(d)), d.iloc[:, 2:].values] @ r[0]
    return np.corrcoef(res.T)[0, 1], len(d)

rows = []
targets = CANDS + ["4EBP1 S65", "PRAS40 S183", "ULK1 S758", "LARP1 S766", "RPS6 S235", "RPS6 S240", "EIF4B S422",
                   "IRS1 S1101", "PRAS40 T246", "GSK3B S9"]
for t in targets:
    y = W[t]
    for ref in R.columns:
        d = pd.concat([y, R[ref]], axis=1).dropna()
        r = d.corr(method="spearman").iloc[0, 1] if len(d) >= 15 else np.nan
        rows.append(dict(site=t, ref=ref, rho=r, n=len(d), kind="spearman"))
    # partials: direct vs S6K vs AKT
    for a, others in [("direct_mTORC1", ["S6K"]), ("S6K", ["direct_mTORC1"]), ("direct_mTORC1", ["AKT"]),
                      ("AKT", ["direct_mTORC1"]), ("ERK_loop", ["PSP_CK2_set"]), ("PSP_CK2_set", ["ERK_loop"]),
                      ("ERK_loop", ["PSP_CDK1_set"]), ("PSP_CDK1_set", ["ERK_loop"])]:
        r, n = pcor(y.rank(), R[a].rank(), [R[o].rank() for o in others])
        rows.append(dict(site=t, ref=f"{a}|{'+'.join(others)}", rho=r, n=n, kind="partial_spearman"))
P = pd.DataFrame(rows)
P.to_csv(f"{OUT}/profile_corr.tsv", sep="\t", index=False)
print(P.pivot_table(index="site", columns="ref", values="rho").round(2).to_string())
print(P.pivot_table(index="site", columns="ref", values="n").to_string())

# ---------- MK2206 timing (A549, 30 vs 60 min, matched doses) ----------
T = L[(L.study_pmid == "36926954") & (L.agents == "MK2206") & L.dose_nM.notna()]
T = T.pivot_table(index=["label", "group", "dose_nM"], columns="time_min", values="log2fc").reset_index()
T.columns = [str(c) for c in T.columns]
T.to_csv(f"{OUT}/mk2206_timing.tsv", sep="\t", index=False)
hi = T[T.dose_nM >= 1000].groupby(["label", "group"])[["30.0", "60.0"]].mean()
hi["frac_at_30"] = hi["30.0"] / hi["60.0"]
print(hi.round(2).to_string())

# ---------- stimulation ----------
St = L[L.drug_class.fillna("").str.startswith("stim") | (L.drug_class == "starvation")]
St = St.groupby(["label", "group", "drug_class", "study_pmid", "cell_line_norm", "time_min"], dropna=False).z.mean().reset_index()
St.to_csv(f"{OUT}/stim_effects.tsv", sep="\t", index=False)
