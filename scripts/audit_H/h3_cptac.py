"""Audit H3: CPTAC fast-vs-stable coupling, adjusted for ischaemia footprint + top-5 phospho PCs, specificity
(coupling to OTHER families' footprints), matching on tumour variance & missingness, cluster bootstrap by site.
Out: results/audit_H/h3_site_cohort.tsv, h3_summary.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import rankdata
sys.path.insert(0, "/workspace/scripts/b_I28"); from lib28 import classify_units, C28, C26
sys.path.insert(0, "/workspace/scripts/b_I26"); from lib26 import gold_sets
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
OUT = "/workspace/results/audit_H"; rng = np.random.default_rng(3)
MAIN = ["AKT", "CDK1/2/5/9", "ERK", "MTOR", "EGFR/ERBB2"]
A = pd.read_csv("/workspace/results/b_I28/tables/turnover_atlas.tsv", sep="\t")
L = pd.read_parquet(f"{C28}/L_units.parquet")
AR = classify_units(L[L.set == "random"])
U = pd.concat([A[A.group.isin(MAIN) & (A.n_units >= 2)][["group", "site", "cls", "R", "gold"]].assign(kind="substrate"),
               AR[AR.group.isin(MAIN) & (AR.n_units >= 2)][["group", "site", "cls", "R"]].assign(kind="random", gold=False)])
E = pd.read_parquet(f"{C26}/psp_edges_annot.parquet"); gold = gold_sets(E)
ISC = pd.read_csv("/workspace/results/val_E/ischaemia_scores.tsv", sep="\t", index_col=0)
ISC = ISC[ISC.variant == "input_tc_padj"].isch


def resid(Y, X):
    """Y: sites x samples (NaN allowed); X: samples x k (no NaN). OLS per site on observed samples (+intercept)."""
    out = np.full(Y.shape, np.nan)
    Xc = np.column_stack([np.ones(len(X)), X])
    for i in range(Y.shape[0]):
        ok = ~np.isnan(Y[i])
        if ok.sum() >= 20:
            b = np.linalg.lstsq(Xc[ok], Y[i, ok], rcond=None)[0]; out[i, ok] = Y[i, ok] - Xc[ok] @ b
    return out


def spearman_rows(Rm, f):
    """Spearman of each row of Rm with vector f (pairwise complete, >=20)."""
    out = np.full(Rm.shape[0], np.nan)
    for i in range(Rm.shape[0]):
        ok = ~np.isnan(Rm[i]) & ~np.isnan(f)
        if ok.sum() >= 20:
            out[i] = np.corrcoef(rankdata(Rm[i, ok]), rankdata(f[ok]))[0, 1]
    return out


rows = []
need = set(U.site) | set().union(*[gold[g] for g in MAIN])
for co in cl.COHORTS:
    t = cl.tumour_ids(co)
    ph = cl.matrix(co, "phospho"); ph = ph[[c for c in ph.columns if c in t]]
    pr = cl.matrix(co, "protein").reindex(columns=ph.columns)
    # global PCs from well-measured phosphosites (>=80% present), z-scored, 0-imputed
    W = ph[ph.notna().mean(axis=1) >= 0.8]; Wz = W.sub(W.mean(1), axis=0).div(W.std(1), axis=0).fillna(0).values.T
    u, s, vt = np.linalg.svd(Wz - Wz.mean(0), full_matrices=False); PC = u[:, :5] * s[:5]
    isc = ISC.reindex(ph.columns)
    if isc.isna().mean() > 0.2:
        print(co, "isch missing frac", isc.isna().mean())
    isc = isc.fillna(isc.median()).values
    P = ph[ph.index.isin(need)]; P = P[P.notna().sum(axis=1) >= 20]
    # protein residualisation (as I28)
    Y = P.values.astype(float); Rp = np.full(Y.shape, np.nan)
    for i, s_ in enumerate(P.index):
        a = s_.split("_")[0]
        if a in pr.index:
            x = pr.loc[a].values.astype(float); ok = ~np.isnan(Y[i]) & ~np.isnan(x)
            if ok.sum() >= 20:
                b = np.polyfit(x[ok], Y[i, ok], 1); Rp[i, ok] = Y[i, ok] - np.polyval(b, x[ok])
    variants = {"I28": Rp, "isch": resid(Rp, isc[:, None]), "isch+5PC": resid(Rp, np.column_stack([isc, PC]))}
    idx = {s_: i for i, s_ in enumerate(P.index)}
    nt = (~np.isnan(Rp)).sum(1)
    for vn, Rm in variants.items():
        Rz = (Rm - np.nanmean(Rm, 1, keepdims=True)) / np.nanstd(Rm, 1, keepdims=True)
        fps = {}
        for g in MAIN:
            gs = [s_ for s_ in gold[g] if s_ in idx and nt[idx[s_]] >= 20]
            if len(gs) < 4: continue
            Gi = [idx[s_] for s_ in gs]; tot = np.nansum(Rz[Gi], 0); cnt = (~np.isnan(Rz[Gi])).sum(0)
            fps[g] = (gs, Gi, tot, cnt)
        for g in fps:
            uu = U[(U.group == g) & U.site.isin(idx)]
            ii = np.array([idx[s_] for s_ in uu.site])
            res = {}
            for g2, (gs, Gi, tot, cnt) in fps.items():
                c = np.full(len(ii), np.nan)
                for k, (s_, i) in enumerate(zip(uu.site, ii)):
                    if s_ in gs:
                        zi = Rz[i]; f = (tot - np.nan_to_num(zi)) / (cnt - (~np.isnan(zi)))
                    else:
                        f = tot / np.maximum(cnt, 1) * np.where(cnt > 0, 1, np.nan)
                    ok = ~np.isnan(Rm[i]) & ~np.isnan(f)
                    if ok.sum() >= 20:
                        c[k] = np.corrcoef(rankdata(Rm[i, ok]), rankdata(f[ok]))[0, 1]
                res[g2] = c
            own = res[g]; others = np.nanmean(np.column_stack([res[g2] for g2 in res if g2 != g]), 1) if len(res) > 1 else np.nan
            rows.append(pd.DataFrame(dict(cohort=co, variant=vn, group=g, site=uu.site.values, kind=uu.kind.values, cls=uu.cls.values,
                                          gold=uu.gold.values, coupling=own, coupling_other=others, n_tum=nt[ii],
                                          sd_resid=np.nanstd(Rp[ii], 1), miss=1 - nt[ii] / Rp.shape[1])))
    print(co, flush=True)
D = pd.concat(rows, ignore_index=True)
D.to_csv(f"{OUT}/h3_site_cohort.tsv", sep="\t", index=False)
