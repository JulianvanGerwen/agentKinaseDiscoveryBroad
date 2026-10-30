"""I03 step 3: tumour co-variation of each site with its candidate kinase's ACTIVATING-site phosphorylation.
Activating sites: PSP Regulatory_sites, human, ON_FUNCTION contains 'enzymatic activity, induced', on the kinase itself.
Data: CPTAC input_tc_padj (site residual after regressing on host protein; tumour-centred), tumours only, 10 cohorts.
Kinase state per tumour = mean of standardised activating-site residuals. Per cohort Pearson r (>=20 tumours);
meta: Stouffer of Fisher z * sqrt(n-3) across cohorts -> tum_z; tum_ncoh = cohorts used; tum_npos = cohorts with r>0.
Outputs cache/b_I03/tumour_long.parquet, cache/b_I03/activating_sites.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
OUT = "/workspace/cache/b_I03"
KC = pd.read_csv(f"{OUT}/kinase_conditions.tsv", sep="\t")
panel = sorted(KC.kinase.unique())

R = pd.read_csv("/workspace/data/Regulatory_sites", sep="\t", skiprows=3, low_memory=False)
R = R[(R.ORGANISM == "human") & R.MOD_RSD.str.endswith("-p", na=False)]
R = R[R.ON_FUNCTION.str.contains("enzymatic activity, induced", na=False)]
R["site"] = R.ACC_ID.str.split("-").str[0] + "_" + R.MOD_RSD.str[:-2]
R = R[~R.ACC_ID.str.contains("-")]
A = R[R.GENE.isin(panel)][["GENE", "site"]].drop_duplicates()
A.to_csv(f"{OUT}/activating_sites.tsv", sep="\t", index=False)
print(A.groupby("GENE").site.apply(list).to_string())

res = {}
for coh in cl.COHORTS:
    X = cl.matrix(coh, "input_tc_padj")
    tum = [c for c in X.columns if c.endswith("_T")]
    X = X[tum]
    X = X[X.notna().sum(1) >= 20]
    Xs = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
    V = Xs.values.astype(np.float64); Mk = ~np.isnan(V)
    for k in panel:
        aset = [s for s in A.site[A.GENE == k] if s in Xs.index]
        if not aset:
            continue
        st = Xs.loc[aset].mean(0, skipna=True).values
        sm = ~np.isnan(st)
        if sm.sum() < 20:
            continue
        M = Mk & sm[None, :]
        n = M.sum(1)
        S = np.where(M, st[None, :], np.nan); Y = np.where(M, V, np.nan)
        mx = np.nanmean(S, 1, keepdims=True); my = np.nanmean(Y, 1, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.nansum((S - mx) * (Y - my), 1) / np.sqrt(np.nansum((S - mx) ** 2, 1) * np.nansum((Y - my) ** 2, 1))
        r = np.where(n >= 20, r, np.nan)
        df = pd.DataFrame(dict(site=Xs.index, r=r, n=n)).dropna()
        df = df[~df.site.isin(aset)]
        res.setdefault(k, []).append(df.assign(cohort=coh))
    print(coh, X.shape)
out = []
for k, L in res.items():
    d = pd.concat(L)
    d["zf"] = np.arctanh(d.r.clip(-0.999, 0.999)) * np.sqrt(d.n - 3)
    g = d.groupby("site").agg(tum_z=("zf", "sum"), tum_ncoh=("zf", "size"), tum_npos=("r", lambda x: (x > 0.1).sum()),
                              tum_rmed=("r", "median"))
    g["tum_z"] = g.tum_z / np.sqrt(g.tum_ncoh)
    out.append(g.reset_index().assign(kinase=k))
T = pd.concat(out, ignore_index=True)
T.to_parquet(f"{OUT}/tumour_long.parquet"); print(T.shape, T.kinase.nunique())
