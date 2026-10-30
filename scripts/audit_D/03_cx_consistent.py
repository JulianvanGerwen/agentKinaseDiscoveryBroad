"""audit_D 03: sites consistently decreased by CX-4945 across Hijazi HL60+NTERA2 (+HeLa), their CK2 motif enrichment,
and specificity vs TBCA / other Hijazi inhibitors; LY2090314 check."""
import common, numpy as np, pandas as pd
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
OUT = "/workspace/results/audit_D/tables"
C = common.conds().set_index("conditionID"); Z = common.zmat()
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
hl, mc, nt, he = "31959955_172_log2", "31959955_176_log2", "31959955_180_log2", "29214152_qPTM_1_log2"
tb = ["31959955_8_log2", "31959955_16_log2", "31959955_24_log2"]
Hcid = C.index[(C.study_pmid == "31959955")]
ZH = Z[Hcid]
# per-site: CX-4945 mean z over HL60/NTERA2 (both measured)
d = Z[[hl, nt, he, mc] + tb].copy()
both = d[[hl, nt]].notna().all(1)
d = d[both]; d["cx_max"] = d[[hl, nt]].max(1)   # both below threshold => max < thr
d["ck2pct"] = M.CSNK2A1.reindex(d.index)
d = d[S.residue.reindex(d.index).isin(["S", "T"])]
for thr in [-1.5, -2, -2.5]:
    sel = d.cx_max < thr
    print(f"both HL60&NTERA2 z<{thr}: n={sel.sum()} of {len(d)}; CK2 motif pct median {d.ck2pct[sel].median():.0f} vs {d.ck2pct.median():.0f};"
          f" frac CK2pct>=90 {(d.ck2pct[sel]>=90).mean():.2f} vs {(d.ck2pct>=90).mean():.2f}")
sel = d.cx_max < -2
L = d[sel].copy(); L["gene"] = S.gene.reindex(L.index); L["window"] = S.sequence_window.reindex(L.index)
L["psp_kin"] = S.psp_kinases_human.reindex(L.index)
# specificity: fraction of other Hijazi inhibitors (same lines) with z<-2
oth = C.index[(C.study_pmid == "31959955") & ~C.agents.isin(["CX4945"])]
L["frac_other_hijazi_lt-2"] = (Z.loc[L.index, oth] < -2).sum(1) / Z.loc[L.index, oth].notna().sum(1)
L = L.sort_values("cx_max"); L.to_csv(f"{OUT}/cx4945_consistent_down_hijazi.tsv", sep="\t")
print(L.round(2).to_string())
# KPNA4 S60 in this framework
k = "O00629_S60"; print("KPNA4 S60 frac other Hijazi inhibitors z<-2:", round(float((Z.loc[k, oth] < -2).sum() / Z.loc[k, oth].notna().sum()), 3))
# LY2090314 targets
Dt = pd.read_csv("/workspace/results/infra_compendium/drug_targets.tsv", sep="\t").set_index("agent")
print(Dt.loc["LY2090314"].to_string())
ly = "31959955_67_log2"; z = Z[ly].dropna()
for nm, mem in [("CK2motif>=95", M.index[M.CSNK2A1 >= 95])]:
    print("LY2090314 HL60", nm, "median z", round(float(z.reindex(mem).dropna().median()), 2), "all", round(float(z.median()), 2))
print("LY2090314 HL60 values for CX-consistent sites:", Z.loc[L.index, ly].round(2).to_dict())
