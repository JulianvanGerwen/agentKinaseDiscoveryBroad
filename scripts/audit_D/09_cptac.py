"""audit_D 09: CPTAC tumour co-variation with a site-level null. For each readout (4EBP1 T37/46, ERK pTpY, CK2 composite, CDK composite),
compute per-cohort Pearson r of every site (protein-adjusted tumour residuals, input_tc_padj), Stouffer across cohorts; report the
candidate's meta-z and its percentile among (a) all sites, (b) motif-matched sites (CK2 pct>=95 / S-P sites / MTOR pct>=75)."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
S = pd.read_csv("/workspace/results/infra_compendium/sites.tsv", sep="\t", low_memory=False).set_index("site")
M = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
pk = S.psp_kinases_human.fillna("")
CK2set = set(S.index[pk.str.contains("CSNK2A1")]) - {"O00629_S60"}
CDKset = set(S.index[pk.str.contains("CDK1")]) - {"Q9Y6I3_S435", "Q9Y6I3_S454"}
READ = {"4EBP1_T37T46": ["Q13541_T37", "Q13541_T46"], "ERK_pTpY": ["P27361_T202", "P27361_Y204", "P28482_T185", "P28482_Y187"]}
CANDS = {"KPNA4 S60": "O00629_S60", "LAMTOR1 S98": "Q6IAA8_S98", "LAMTOR1 S97": "Q6IAA8_S97", "RRAGC S15": "Q9HB90_S15",
         "EPN1 S435": "Q9Y6I3_S435", "PRAS40 S183": "Q96B36_S183", "4EBP1 S65": "Q13541_S65", "LARP1 S766": "Q6PKG0_S766",
         "RSK1 S380": "Q15418_S380", "NPM1 S125": "P06748_S125", "EIF2S2 S2": "P20042_S2", "SF3A1 S451": "Q15459_S451", "ANLN S54": "Q9NQW6_S54"}
acc = {}; cnt = {}
for coh in cl.COHORTS:
    X = cl.matrix(coh, "input_tc_padj"); X = X[[c for c in X.columns if c.endswith("_T")]]
    X = X[X.notna().sum(1) >= 30]
    Xs = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
    R = {k: Xs.loc[[s for s in v if s in Xs.index]].mean(0) for k, v in READ.items() if any(s in Xs.index for s in v)}
    for nm, st in [("CK2_comp", CK2set), ("CDK_comp", CDKset)]:
        ss = [s for s in st if s in Xs.index and Xs.loc[s].notna().mean() > .5]
        if len(ss) >= 5: R[nm] = Xs.loc[ss].mean(0)
    for k, r in R.items():
        A = Xs.values; m = ~np.isnan(A) & ~np.isnan(r.values)[None, :]
        n = m.sum(1); a = np.where(m, A, 0); b = np.where(m, r.values[None, :], 0)
        ma = a.sum(1) / n; mb = b.sum(1) / n
        cov = (np.where(m, a - ma[:, None], 0) * np.where(m, b - mb[:, None], 0)).sum(1)
        va = (np.where(m, a - ma[:, None], 0) ** 2).sum(1); vb = (np.where(m, b - mb[:, None], 0) ** 2).sum(1)
        rr = cov / np.sqrt(va * vb); zf = np.arctanh(np.clip(rr, -.999, .999)) * np.sqrt(np.maximum(n - 3, 1))
        zf = pd.Series(np.where(n >= 20, zf, np.nan), index=Xs.index)
        for key in (k,):
            acc.setdefault(key, []).append(zf)
    print(coh, list(R))
rows = []
for k, lst in acc.items():
    D = pd.concat(lst, axis=1); nco = D.notna().sum(1); meta = D.sum(1) / np.sqrt(nco); meta = meta[nco >= 3]
    for lab, s in CANDS.items():
        if s not in meta.index: continue
        r = dict(site=lab, readout=k, meta_z=meta[s], n_coh=int(nco[s]), pct_all=(meta < meta[s]).mean() * 100)
        for nm, mem in [("pct_CK2motif95", M.index[M.CSNK2A1 >= 95]), ("pct_SP", S.index[S.sequence_window.str[8] == "P"]),
                        ("pct_MTOR75", M.index[M.MTOR >= 75])]:
            mm = meta.reindex(mem).dropna(); r[nm] = (mm < meta[s]).mean() * 100
        rows.append(r)
    print(k, "null meta-z median %.2f, 95th %.2f" % (meta.median(), meta.quantile(.95)))
T = pd.DataFrame(rows); T.to_csv("/workspace/results/audit_D/tables/cptac_null.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(T.round(1).to_string(index=False))
