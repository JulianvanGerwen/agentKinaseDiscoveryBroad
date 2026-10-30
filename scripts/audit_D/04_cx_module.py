"""audit_D 04: is the CX-4945-specific decreased set (HL60&NTERA2 z<-2, rarely down with other Hijazi inhibitors) enriched for
CK2 motif (vs off-target CLK/DYRK motifs)? Uses table from 03."""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
import common
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
Z = common.zmat(); C = common.conds().set_index("conditionID")
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
hl, nt, he = "31959955_172_log2", "31959955_180_log2", "29214152_qPTM_1_log2"
oth = C.index[(C.study_pmid == "31959955") & ~C.agents.isin(["CX4945"])]
d = Z[[hl, nt, he]].copy(); d = d[d[[hl, nt]].notna().all(axis=1) & S.residue.reindex(d.index).isin(["S", "T"]).values]
d["cx_max"] = d[[hl, nt]].max(axis=1)
d["frac_other"] = (Z.loc[d.index, oth] < -2).sum(axis=1) / Z.loc[d.index, oth].notna().sum(axis=1)
for k in ["CSNK2A1", "CLK1", "CLK3", "DYRK1A", "SRPK1", "HIPK2", "CDK1", "PLK1", "CSNK1D"]:
    d[k] = M[k].reindex(d.index)
spec = (d.cx_max < -2) & (d.frac_other < 0.06)
nons = (d.cx_max < -2) & (d.frac_other >= 0.06)
print("CX-specific down n=%d, CX-down-but-promiscuous n=%d, background n=%d" % (spec.sum(), nons.sum(), len(d)))
for k in ["CSNK2A1", "CLK1", "CLK3", "DYRK1A", "SRPK1", "HIPK2", "CDK1", "PLK1", "CSNK1D"]:
    a = (d.loc[spec, k] >= 90).sum(); b = spec.sum() - a; c_ = (d.loc[~spec, k] >= 90).sum(); e = (~spec).sum() - c_
    orr, p = fisher_exact([[a, b], [c_, e]])
    print(f"{k:8s} spec frac>=90 {a/spec.sum():.2f} | promiscuous {(d.loc[nons,k]>=90).mean():.2f} | bg {(d[k]>=90).mean():.2f} | OR {orr:.1f} p={p:.2g}")
L = d[spec].sort_values("cx_max"); L["gene"] = S.gene.reindex(L.index); L["window"] = S.sequence_window.reindex(L.index)
print(L[["gene", "window", hl, nt, he, "frac_other", "CSNK2A1", "CLK1"]].round(2).to_string())
# HeLa (verified CK2 inhibition): among the CX-specific Hijazi set with HeLa data, how many fall?
h = L[he].dropna(); print("HeLa z for Hijazi-CX-specific sites (n=%d): median %.2f, frac<-2 %.2f; background frac<-2 %.2f" % (len(h), h.median(), (h < -2).mean(), (Z[he].dropna() < -2).mean()))
L.to_csv("/workspace/results/audit_D/tables/cx4945_specific_module.tsv", sep="\t")
