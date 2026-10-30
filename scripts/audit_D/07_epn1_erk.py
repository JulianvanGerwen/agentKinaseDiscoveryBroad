"""audit_D 07: EPN1 S435 vs ERK pathway. Unit-level concordance with ERK loop (MEKi/RAFi/EGFRi) or RSK1 S363 (ERKi),
LOSO, null of matched random S/T-P sites, same-protein controls; KPNA4 S60 in the same ERK-verified units."""
import common, numpy as np, pandas as pd
from scipy.stats import spearmanr
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
OUT = "/workspace/results/audit_D/tables"
U = pd.read_parquet("/workspace/cache/audit_D/unit_z.parquet")
ui = pd.read_csv("/workspace/cache/audit_D/unit_info.tsv", sep="\t", index_col=0, dtype={"study_pmid": str})
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
sel = ui.index[ui.drug_class.isin(["MEKi", "ERKi", "RAFi", "EGFRi"])]
Us = U[sel]
loop = Us.loc[[s for s in ["P27361_T202", "P27361_Y204", "P28482_T185", "P28482_Y187"] if s in Us.index]].mean()
rsk = Us.loc["Q15418_S363"]
down = pd.Series(np.where(ui.loc[sel, "drug_class"] == "ERKi", rsk < -2, loop < -2), index=sel)
known = pd.Series(np.where(ui.loc[sel, "drug_class"] == "ERKi", rsk.notna(), loop.notna()), index=sel)
print("units:", len(sel), "pathway verifiable:", known.sum(), "verified down:", (down & known).sum())
def stats(s):
    v = Us.loc[s]; a = v[down & known].dropna(); b = v[~down & known].dropna()
    d = pd.concat([loop, v], axis=1).dropna()
    rho = spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0] if len(d) > 8 else np.nan
    return dict(n_down_units=len(a), med_down=a.median(), n_lt2_down=int((a < -2).sum()), n_notdown=len(b), med_notdown=b.median(),
                rho_loop=rho, n_rho=len(d), n_studies=ui.loc[a.index, "study_pmid"].nunique())
rows = {}
for s in ["Q9Y6I3_S435", "Q9Y6I3_S454", "Q9Y6I3_T470", "Q9Y6I3_S447", "Q9Y6I3_T460", "O00629_S60", "Q15418_S380", "Q15418_S363", "Q6IAA8_S98"]:
    if s in Us.index: rows[S.loc[s, "gene"] + " " + s.split("_")[1]] = stats(s)
R = pd.DataFrame(rows).T; print(R.round(2).to_string())
# LOSO for EPN1 S435
v = Us.loc["Q9Y6I3_S435"]; st = ui.loc[sel, "study_pmid"]
lo = []
for x in st[down & known & v.notna()].unique():
    m = (st != x); a = v[down & known & m].dropna(); d = pd.concat([loop[m], v[m]], axis=1).dropna()
    lo.append(dict(left_out=x, n=len(a), med_down=a.median(), rho=spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0]))
L = pd.DataFrame(lo); print(L.round(2).to_string())
# null: S/T-P sites measured in >= 70% as many verified-down units
stp = S.index[(S.sequence_window.str[8] == "P") & S.residue.isin(["S", "T"])].intersection(Us.index)
X = Us.loc[stp, (down & known)[lambda x: x].index]
ok = X.notna().sum(axis=1) >= 0.7 * R.loc["EPN1 S435", "n_down_units"]
med = X[ok].median(axis=1)
print("null S/T-P sites n=%d: EPN1 S435 median-in-down-units %.2f, percentile %.1f; frac null <= -1: %.3f" % (
    ok.sum(), R.loc["EPN1 S435", "med_down"], (med < R.loc["EPN1 S435", "med_down"]).mean() * 100, (med <= -1).mean()))
# null rho
rhos = []
for s in X[ok].index:
    d = pd.concat([loop, Us.loc[s]], axis=1).dropna()
    if len(d) >= 20: rhos.append(spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0])
rhos = np.array(rhos); print("null rho (S/T-P) n=%d: median %.2f, 95th %.2f, 99th %.2f; EPN1 S435 rho %.2f -> pct %.1f" % (
    len(rhos), np.median(rhos), np.quantile(rhos, .95), np.quantile(rhos, .99), R.loc["EPN1 S435", "rho_loop"], (rhos < R.loc["EPN1 S435", "rho_loop"]).mean() * 100))
# per unit list
T = pd.DataFrame({"ERK_loop": loop, "RSK1_S363": rsk, "EPN1_S435": Us.loc["Q9Y6I3_S435"], "EPN1_S454": Us.loc["Q9Y6I3_S454"],
                  "KPNA4_S60": Us.loc["O00629_S60"], "down": down, "known": known}).join(ui).sort_values("ERK_loop")
T.to_csv(f"{OUT}/epn1_erk_units.tsv", sep="\t"); print(T[T.EPN1_S435.notna() | T.KPNA4_S60.notna()].round(2).to_string())
