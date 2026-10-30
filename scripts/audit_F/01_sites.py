"""Task 1a/1b: per-site, per-inhibitor-class values with own canonical-marker verification; lab-level counts; leave-one-lab-out."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__)); from lib_f import *

OUT = "/workspace/results/audit_F/tables"; os.makedirs(OUT, exist_ok=True)
SITES = {  # site: (gene, psp classes, claimed alt)
    "O14757_S280": ("CHEK1", ["AKT"], "ERK"), "P29317_S897": ("EPHA2", ["AKT"], "MEK"), "Q14247_S405": ("CTTN", ["PKC", "CDK"], "MEK"),
    "P36507_T394": ("MAP2K2", ["CDK"], "RAF"), "P26358_S143": ("DNMT1", ["AKT"], "MEK"), "O00443_S259": ("PIK3C2A", ["CDK"], "MEK"),
    "Q01167_S373": ("FOXK2", ["CDK"], "MTOR"), "Q96F86_S161": ("EDC3", ["AKT"], "ALK"), "Q9UJU6_S269": ("DBNL", ["AKT"], "MEK"),
    "P0C1Z6_S180": ("TFPT", ["CDK"], "RAF"), "Q9UKA4_T1100": ("AKAP11", ["MTOR"], "ERK"),
}
CC = conditions()
Z = zmatrix(CC.conditionID.unique())
ms = marker_score(Z, CC)
CC["marker_z"] = [ms[(c, i)][0] for c, i in zip(CC.cls, CC.conditionID)]
CC["marker_n"] = [ms[(c, i)][1] for c, i in zip(CC.cls, CC.conditionID)]
CC["verified"] = CC.marker_z <= -2.5
# cross-talk readouts in every condition
sgn = np.where(CC.rev, -1, 1)
for m in ["AKT", "MAPK", "MTOR"]:
    v = Z.reindex(MARK[m]).T
    CC[f"{m}_mark"] = [(v.loc[c].dropna().sum() / np.sqrt(max(v.loc[c].notna().sum(), 1))) * s if v.loc[c].notna().sum() >= 2 else np.nan
                       for c, s in zip(CC.conditionID, sgn)]
CC.to_csv(f"{OUT}/conditions_classes.tsv", sep="\t", index=False, float_format="%.3g")
print(CC.groupby("cls").agg(n=("verified", "size"), ver=("verified", "sum"), labs=("lab", "nunique"),
                            labs_ver=("lab", lambda s: s[CC.loc[s.index, "verified"]].nunique())).to_string())

TEST = ["AKT", "PKC", "CDK", "MTOR", "MEK", "ERK", "RAF", "RSK", "ALK", "EGFR"]
rows = []
for s, (g, psp, alt) in SITES.items():
    for r in CC[CC.cls.isin(TEST)].itertuples():
        if s not in Z.index or np.isnan(Z.at[s, r.conditionID]):
            continue
        rows.append(dict(site=s, gene=g, role="psp" if r.cls in psp else ("claimed_alt" if r.cls == alt else "other"), cls=r.cls,
                         conditionID=r.conditionID, study=r.study, lab=r.lab, agent=r.agent, cell=r.cell, time_min=r.time_min,
                         z=float(Z.at[s, r.conditionID]) * (-1 if r.rev else 1), marker_z=r.marker_z, verified=r.verified,
                         AKT_mark=r.AKT_mark, MAPK_mark=r.MAPK_mark, MTOR_mark=r.MTOR_mark))
Lg = pd.DataFrame(rows)
Lg.to_csv(f"{OUT}/site_condition_long.tsv", sep="\t", index=False, float_format="%.3g")

# lab-level summary under verified inhibition
MAPKC = ["MEK", "ERK", "RAF", "RSK"]


def lab_summary(d):
    d = d[d.verified]
    if d.empty:
        return dict(n_lab=0, lab_dec=0, lab_flat=0, labs="", lab_z="")
    st = d.groupby(["lab", "study"]).z.mean().groupby("lab").mean()
    return dict(n_lab=len(st), lab_dec=int((st <= -1).sum()), lab_flat=int((st > -0.5).sum()), labs=";".join(st.index),
                lab_z=";".join(f"{l}:{v:.2f}" for l, v in st.items()), n_cond=len(d))


summ = []
for s, (g, psp, alt) in SITES.items():
    d = Lg[Lg.site == s]
    for lab_name, dd in [("PSP:" + "/".join(psp), d[d.cls.isin(psp)]), ("MAPK(MEK/ERK/RAF/RSK)", d[d.cls.isin(MAPKC)]),
                         ("claimed_alt:" + alt, d[d.cls == alt])] + [(c, d[d.cls == c]) for c in ["RSK", "ALK", "EGFR", "MTOR"] if c != alt]:
        x = lab_summary(dd)
        # unverified, for transparency
        uv = dd[~dd.verified]
        x.update(site=s, gene=g, test=lab_name, n_cond_unverified=len(uv),
                 unverified_lab_z=";".join(f"{l}:{v:.2f}" for l, v in uv.groupby("lab").z.mean().items()))
        # leave-one-lab-out: does 'lab_dec>=1 and majority dec' hold without each lab?
        ver = dd[dd.verified]
        labs = ver.lab.unique()
        lolo = []
        for L in labs:
            st = ver[ver.lab != L].groupby("lab").z.mean()
            lolo.append(f"-{L}:{(st <= -1).sum()}/{len(st)}")
        x["lolo_dec_over_labs"] = ";".join(lolo)
        summ.append(x)
S = pd.DataFrame(summ)[["gene", "site", "test", "n_lab", "lab_dec", "lab_flat", "lab_z", "n_cond", "lolo_dec_over_labs", "n_cond_unverified", "unverified_lab_z"]]
S.to_csv(f"{OUT}/site_lab_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(S.to_string())
