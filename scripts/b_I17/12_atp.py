"""ATP-competition (Cheng-Prusoff) test of the potency shift (PI follow-up).
Per drug: median dpEC50 of TARGET-SUBSTRATE sites only (down-regulated; on target-cluster proteins, curated
PSP/OmniPath substrate of T, or kinaseSub top5 predicted for T); median over experiments of per-experiment medians.
Annotations: Klaeger S1 binding type (1 type I, 2 type II, 3 allosteric type III, 4 allosteric type IV, 0 unknown),
binding mode (covalent); native ATP Kd,app (uM) from Patricelli et al. 2011 Chem Biol Table 1 (KiNativ, PMID 21700206);
Cheng-Prusoff predicted shift = -log10(1+[ATP]/Km) for [ATP] = 1 and 5 mM.
Second affinity reference: Davis et al. 2011 KINOMEscan Kd (ATP-free, recombinant) where available.
Tests: (a) ATP-site (types 1/2) vs non-ATP (3/4) Mann-Whitney + leave-one-drug-out; (b) Spearman shift vs Km;
(c) type I vs type II. Output: results/b_I17/atp_drug_table.tsv, atp_tests.tsv, fig6_atp.png
"""
import numpy as np, pandas as pd, glob, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
RES = "/workspace/results/b_I17"
A = pd.read_csv(f"{RES}/shift_sites.tsv", sep="\t")
A["tsub"] = A.self | A.ks_direct | A.pred_direct
dn = A[(A.direction == "down")]
s1 = pd.read_excel("/workspace/data/Klaeger/aan4368_Table_S1.xlsx", sheet_name="Inhibitor annotation").set_index("Drug")
DM = pd.read_csv(f"{RES}/drug_map.tsv", sep="\t").set_index("agent")
# Patricelli 2011 Table 1 native ATP Kd,app (uM); others not in table -> NaN (not guessed)
KM = {"EGFR": 18.0, "SRC": 9.8, "ABL1": 5.2, "MAP2K1": 8.1}
# Davis 2011 Kd (nM), non-mutant entries
dv = pd.read_excel(glob.glob("<HACKATHON_ROOT>/kinaseSubMain/external/davis2011/*MOESM5*")[0])
dv = dv[~dv.Kinase.astype(str).str.contains(r"\(")]
dcols = {c.lower().replace("-", ""): c for c in dv.columns}
rows = []
for (ex, ag), d in dn.groupby(["exp_id", "agent"]):
    t = d[d.tsub]
    rows.append(dict(exp_id=ex, agent=ag, target=d.target.iloc[0], pKd=d.pKd_target.iloc[0], n_tsub=len(t),
                     shift_tsub=t["shift"].median() if len(t) else np.nan, shift_all=d["shift"].median(),
                     pEC50_tsub=t.pEC50.median() if len(t) else np.nan))
E = pd.DataFrame(rows)
D = E.groupby("agent").agg(target=("target", "first"), pKd=("pKd", "first"), n_exp=("exp_id", "size"), n_tsub=("n_tsub", "sum"),
                           shift_tsub=("shift_tsub", "median"), shift_all=("shift_all", "median"), pEC50_tsub=("pEC50_tsub", "median")).reset_index()
D["klaeger"] = D.agent.map(DM.klaeger)
D["binding_type"] = [int(s1.loc[k, "Binding type"]) for k in D.klaeger]
D["covalent"] = [s1.loc[k, "Binding mode"] == "covalent" for k in D.klaeger]
D["atp_site"] = D.binding_type.isin([1, 2]) | ((D.binding_type == 0) & (D.agent == "CUDC101"))  # CUDC-101 is an ATP-site quinazoline
D["Km_ATP_uM"] = D.target.map(KM)
for atp in [1000, 5000]:
    D[f"CP_pred_{atp//1000}mM"] = np.where(D.atp_site, -np.log10(1 + atp / D.Km_ATP_uM), 0.0)
DAVIS_ALIAS = {"selumetinib": "AZD-6244/ARRY-886", "afatinib": "BIBW-2992", "pictilisib": "GDC-0941"}  # BIBF-1120 'derivative' not used for nintedanib
def davis_kd(agent, target):
    k = DAVIS_ALIAS.get(agent) or dcols.get(agent.lower().replace("-", ""))
    if k is None: return np.nan
    v = dv.loc[dv["Entrez Gene Symbol"] == target, k]
    v = pd.to_numeric(v, errors="coerce").dropna()
    return v.min() if len(v) else np.nan
D["davis_kd_nM"] = [davis_kd(a, t) for a, t in zip(D.agent, D.target)]
D["shift_tsub_davis"] = D.pEC50_tsub - (9 - np.log10(D.davis_kd_nM))
D.to_csv(f"{RES}/atp_drug_table.tsv", sep="\t", index=False)
print(D.round(2).to_string())

T = []
def mwu(col, sub=None):
    d = D if sub is None else D[sub]
    d = d.dropna(subset=[col])
    a, b = d.loc[~d.atp_site, col], d.loc[d.atp_site, col]
    if len(a) < 2 or len(b) < 2: return None
    p = mannwhitneyu(a, b, alternative="greater").pvalue
    lodo = []
    for ag in d.agent:
        dd = d[d.agent != ag]
        lodo.append(dd.loc[~dd.atp_site, col].median() - dd.loc[dd.atp_site, col].median())
    return dict(test=f"(a) non-ATP-site minus ATP-site, {col}", n_nonATP=len(a), n_ATP=len(b), median_nonATP=a.median(),
                median_ATP=b.median(), effect=a.median() - b.median(), p=p, lodo_min=min(lodo), lodo_max=max(lodo))
for col in ["shift_tsub", "shift_all", "shift_tsub_davis"]:
    r = mwu(col)
    if r: T.append(r)
r = mwu("shift_tsub", ~D.agent.isin(["PD0325901", "refametinib", "selumetinib"]))  # MK2206-only allosteric arm
if r: r["test"] += " [excluding MEKi]"; T.append(r)
k = D[D.atp_site & D.Km_ATP_uM.notna() & ~D.covalent].dropna(subset=["shift_tsub"])
if len(k) >= 3:
    rho, p = spearmanr(k.Km_ATP_uM, k.shift_tsub)
    T.append(dict(test="(b) Spearman shift_tsub vs Km(ATP), ATP-site non-covalent", n_ATP=len(k), effect=rho, p=p,
                  note=";".join(f"{a}:{km}" for a, km in zip(k.agent, k.Km_ATP_uM))))
    T.append(dict(test="(b) median observed minus Cheng-Prusoff(1 mM) predicted", n_ATP=len(k),
                  effect=(k.shift_tsub - k.CP_pred_1mM).median(), note="positive = smaller shift than predicted"))
    T.append(dict(test="(b) implied [ATP]/Km from observed shift (median)", n_ATP=len(k), effect=np.median(10 ** (-k.shift_tsub) - 1)))
t1 = D[(D.binding_type == 1)].shift_tsub.dropna(); t2 = D[(D.binding_type == 2)].shift_tsub.dropna()
T.append(dict(test="(c) type II minus type I, shift_tsub", n_nonATP=len(t2), n_ATP=len(t1), effect=t2.median() - t1.median(),
              p=mannwhitneyu(t2, t1).pvalue if len(t2) >= 2 and len(t1) >= 2 else np.nan))
TT = pd.DataFrame(T); TT.to_csv(f"{RES}/atp_tests.tsv", sep="\t", index=False); print(TT.round(3).to_string())

fig, ax = plt.subplots(figsize=(6.4, 4.2))
col = {1: "#2a78d6", 2: "#1baf7a", 3: "#eb6834", 4: "#e87ba4", 0: "#8a8984"}
lab = {1: "type I", 2: "type II", 3: "allosteric (MEK, type III)", 4: "allosteric (type IV)", 0: "ATP-site, unannotated"}
d = D.dropna(subset=["shift_tsub"]).sort_values("shift_tsub").reset_index(drop=True)
for i, r in d.iterrows():
    ax.plot(r.shift_tsub, i, "o", color=col[r.binding_type], ms=7, mec="white", mew=1.2, label=lab[r.binding_type])
    if np.isfinite(r.CP_pred_1mM) and r.atp_site:
        ax.plot([r.CP_pred_5mM, r.CP_pred_1mM], [i, i], color="#8a8984", lw=3, alpha=.5, solid_capstyle="round")
    if np.isfinite(r.shift_tsub_davis):
        ax.plot(r.shift_tsub_davis, i, "D", color="white", mec=col[r.binding_type], mew=1.5, ms=6, label="vs Davis 2011 Kd (ATP-free)")
ax.set_yticks(range(len(d))); ax.set_yticklabels([f"{a}{' (cov.)' if c else ''} → {t}" for a, c, t in zip(d.agent, d.covalent, d.target)], fontsize=7)
ax.axvline(0, color="#52514e", lw=.8, ls="--"); ax.set_xlabel("ΔpEC50 of target-substrate sites (● Kinobeads Kd, ◇ Davis Kd)")
h, l = ax.get_legend_handles_labels(); u = dict(zip(l, h)); ax.legend(u.values(), u.keys(), fontsize=6.5, frameon=False, loc="lower right")
ax.set_title("Grey bar = Cheng–Prusoff shift at 1–5 mM ATP\n(native ATP Kd,app from Patricelli 2011)", fontsize=7.5, loc="left")
fig.tight_layout(); fig.savefig(f"{RES}/fig6_atp.png", dpi=150)
