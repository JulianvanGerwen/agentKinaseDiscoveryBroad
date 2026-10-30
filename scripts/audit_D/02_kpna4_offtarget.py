"""audit_D 02: which kinase motif explains sites decreased by CX-4945 (AUROC of motif pct for z<-2 vs rest);
KPNA4 S60 across all Hijazi inhibitors (specificity of CX-4945 within study); inhibitors hitting CX-4945 off-targets;
KPNA4 S60 dynamics across the compendium."""
import common, numpy as np, pandas as pd
from scipy.stats import rankdata
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
OUT = "/workspace/results/audit_D/tables"
C = common.conds().set_index("conditionID"); Z = common.zmat()
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
K = "O00629_S60"
st = S.residue.isin(["S", "T"])
cx = ["31959955_172_log2", "31959955_176_log2", "31959955_180_log2", "29214152_qPTM_1_log2", "31959955_8_log2", "31959955_16_log2", "31959955_24_log2"]
def auroc(score, lab):
    r = rankdata(score); n1 = lab.sum(); n0 = len(lab) - n1
    return (r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
rows = []
for c in cx:
    z = Z[c].dropna(); z = z[z.index.isin(S.index[st])]; z = z[z.index.isin(M.index)]
    lab = (z < -2).values; mm = M.loc[z.index]
    for k in M.columns:
        rows.append(dict(cid=c, kinase=k, auroc=auroc(mm[k].values.astype(float), lab), n_down=int(lab.sum())))
A = pd.DataFrame(rows).pivot_table(index="kinase", columns="cid", values="auroc")
A = A[[c for c in cx]]
A["mean_CX_hijazi"] = A[cx[:3]].mean(1)
A.to_csv(f"{OUT}/motif_auroc_decreased_ck2i.tsv", sep="\t")
print("n_down:", pd.DataFrame(rows).groupby("cid").n_down.first().to_dict())
for col in ["29214152_qPTM_1_log2", "mean_CX_hijazi"]:
    print("\nTop motifs (AUROC of motif pct for z<-2) sorted by", col); print(A.sort_values(col, ascending=False).head(15).round(3).to_string())
print(A.loc[["CSNK2A1", "CSNK2A2", "CLK1", "CLK3", "DYRK1A", "MAPKAPK3", "PIM1", "GSK3B", "HIPK2", "PLK1", "CDK5", "SLK", "STK10", "BMPR2", "GRK2", "CSNK1D"]].round(3).to_string())

# KPNA4 S60 in every condition where measured
k = Z.loc[K].dropna()
D = C.loc[k.index, ["study_pmid", "agents", "cell_line_norm", "dose_nM", "time_min", "perturbation_class", "inhibitor_profiled_targets"]].copy()
D["z"] = k; D["lfc"] = common.lfc(sites=[K]).loc[K, k.index]
D = D.sort_values("z"); D.to_csv(f"{OUT}/kpna4_s60_all_conditions.tsv", sep="\t")
print("\nKPNA4 S60 measured in", len(D), "conditions,", D.study_pmid.nunique(), "studies; |z|>2 frac", round((D.z.abs() > 2).mean(), 3))
print(D.head(30).drop(columns="inhibitor_profiled_targets").round(2).to_string())
# within Hijazi: rank of CX-4945 among inhibitors per line
H = D[D.study_pmid == "31959955"]
for cl, g in H.groupby("cell_line_norm"):
    g = g.sort_values("z"); print(cl, "n inhibitors", len(g), "CX4945 rank", list(g.agents).index("CX4945") + 1 if "CX4945" in list(g.agents) else None,
                                 "| lowest 5:", list(zip(g.agents[:5], g.z.round(2)[:5])))
# inhibitors that hit CX-4945 off-targets (Klaeger <=1uM / Hijazi target lists), excluding CK2
Dt = pd.read_csv("/workspace/results/infra_compendium/drug_targets.tsv", sep="\t")
off = ["CLK1", "CLK2", "CLK3", "CLK4", "DYRK1A", "DYRK1B", "DYRK2", "HIPK1", "HIPK2", "HIPK3", "PIM1", "PIM3", "MAPKAPK3", "PLK1", "GSK3A", "GSK3B", "BMPR2"]
rows = []
for _, r in Dt.iterrows():
    tg = set(str(r.klaeger_kd_le1000nM).replace(",", ";").split(";")) | set(str(r.hijazi_targets).split(";")) | set(str(r.primary_targets).split(";"))
    hit = sorted(tg & set(off)); ck = bool(tg & {"CSNK2A1", "CSNK2A2"})
    if hit: rows.append(dict(agent=r.agent, offt_hits=";".join(hit), hits_CK2=ck))
O = pd.DataFrame(rows)
J = D.reset_index().merge(O, left_on="agents", right_on="agent")
J = J[~(J.dose_nM < 100)]
print("\nKPNA4 S60 under inhibitors hitting CX-4945 off-targets:"); print(J[["index", "study_pmid", "agents", "cell_line_norm", "dose_nM", "z", "offt_hits", "hits_CK2"]].sort_values("z").round(2).to_string())
J.to_csv(f"{OUT}/kpna4_offtarget_inhibitors.tsv", sep="\t", index=False)
# dynamics: fraction of conditions with |z|>2, KPNA4 S60 vs PSP CK2 sites and all sites with >=50 measurements
nm = Z.notna().sum(1); fr = (Z.abs() > 2).sum(1) / nm
pk = S.psp_kinases_human.fillna("")
ck2 = S.index[pk.str.contains("CSNK2A")]
el = nm[nm >= 50].index
print("\nresponsiveness frac|z|>2: KPNA4 S60 %.3f (n=%d); PSP CK2 sites median %.3f (n=%d); all sites median %.3f; KPNA4 S60 pct among all %.1f; among CK2 %.1f" % (
    fr[K], nm[K], fr.reindex(ck2.intersection(el)).median(), len(ck2.intersection(el)), fr[el].median(), (fr[el] < fr[K]).mean() * 100,
    (fr.reindex(ck2.intersection(el)) < fr[K]).mean() * 100))
# CK2-motif (>=95) acidic sites responsiveness vs others
ck2m = M.index[M.CSNK2A1 >= 95].intersection(el)
print("CK2-motif>=95 sites median frac|z|>2: %.3f (n=%d) vs non: %.3f" % (fr[ck2m].median(), len(ck2m), fr[el.difference(ck2m)].median()))
