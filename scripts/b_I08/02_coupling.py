#!/usr/bin/env python3
"""I08 step 2: per cohort x kinase coupling of substrate-inferred activity (kinase-sites masked) with the kinase's
own mRNA, protein and own phosphosites (raw = tc log2; padj = residual on own protein). Spearman (rank) r, and
partial r adjusting purity (WES), ESTIMATE score and proliferation (E2F+G2M ssGSEA).
Nulls: (i) other kinases' protein -> percentile of own protein among all kinase proteins;
       (ii) meanz only: 200 random same-size/same-residue substrate sets -> null mean/sd of r.
Output: results/b_I08/cache/coupling_long.parquet, joint_r2.parquet"""
import os, sys, numpy as np, pandas as pd
os.environ['OMP_NUM_THREADS'] = '8'
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *

kin = pd.read_csv(f'{R}/model_corev3/kinases.tsv', sep='\t')
k2acc = dict(zip(kin.kinase, kin.uniprot))
ann = pd.read_csv(f'{OUT}/tables/kinase_site_annotation.tsv', sep='\t')
# canonical activating site per kinase: PSP-activating or A-loop, most cohorts; prefer PSP-activating & A-loop
cand = ann[ann.psp_activating | ann.in_aloop].copy()
cand['score'] = cand.n_cohorts * 10 + cand.psp_activating * 2 + cand.in_aloop * 1
canon = cand.sort_values('score', ascending=False).drop_duplicates('kinase').set_index('kinase').site_key.to_dict()
pd.Series(canon, name='canonical_site').rename_axis('kinase').to_csv(f'{OUT}/tables/canonical_site.tsv', sep='\t')
SOURCES = ['meanz_tc_padj', 'corev3_tc_padj', 'meanz_tc', 'corev3_tc']


def rk_rows(M):
    return np.vstack([rank(r) for r in M]) if len(M) else M


rows, jrows = [], []
for src in SOURCES:
    meth, var = src.split('_', 1)
    ACT = pd.read_parquet(f'{OUT}/cache/act_{meth}_{var}.parquet')
    for c in cl.COHORTS:
        T = [t for t in cl.tumour_ids(c) if t in ACT.columns]
        C = covariates(c, T); Cr = np.column_stack([rank(C[x]) for x in C])
        prot = cl.matrix(c, 'protein').reindex(columns=T)
        rna = cl.matrix(c, 'rna'); rna = rna[~rna.index.duplicated()].reindex(columns=T)
        raw = cl.matrix(c, 'input_tc').reindex(columns=T); padj = cl.matrix(c, 'input_tc_padj').reindex(columns=T)
        kin_m = kin[kin.uniprot.isin(prot.index)]
        PK = prot.reindex(kin_m.uniprot).values; PKr = rk_rows(PK)
        PKr_p = np.vstack([resid(r, Cr) for r in PKr])
        null = None
        if meth == 'meanz':
            nz = np.load(f'{OUT}/cache/null_{var}_{c}.npz', allow_pickle=True)
            ncols = list(nz['samples']); sel = [ncols.index(t) for t in T]
            null = {k: nz['null'][i][:, sel] for i, k in enumerate(nz['kin'])}
        A = ACT.reindex(columns=T)
        site_by_acc = {}
        for s in raw.index:
            site_by_acc.setdefault(s.split('_')[0], []).append(s)
        for k in A.index:
            if k not in k2acc: continue
            a = A.loc[k].values.astype(float)
            if np.isfinite(a).sum() < 20 or np.nanstd(a) == 0: continue
            acc = k2acc[k]; ar = rank(a); ap = resid(ar, Cr)
            feats = {}
            if acc in prot.index: feats['prot'] = prot.loc[acc].values.astype(float)
            if k in rna.index: feats['rna'] = rna.loc[k].values.astype(float)
            for s in site_by_acc.get(acc, []):
                feats[f'{s}|raw'] = raw.loc[s].values.astype(float)
                if s in padj.index: feats[f'{s}|padj'] = padj.loc[s].values.astype(float)
            if null is not None and k in null:
                N = null[k].astype(float); N = np.where(np.isfinite(N), N, np.nanmean(N, 1, keepdims=True))
                Nr = rk_rows(N); Np = np.vstack([resid(r, Cr) for r in Nr])
            else: Nr = None
            for f, x in feats.items():
                if np.isfinite(x).sum() < 20: continue
                xr = rank(x); xp = resid(xr, Cr)
                r, n = pcor(ar, xr); rp, _ = pcor(ap, xp)
                rec = dict(source=src, cohort=c, kinase=k, feature=f, n=n, r=r, r_partial=rp)
                if Nr is not None:
                    m = np.isfinite(xr)
                    nr = np.array([np.corrcoef(Nr[i, m], xr[m])[0, 1] for i in range(len(Nr))])
                    mp = np.isfinite(xp)
                    npr = np.array([np.corrcoef(Np[i, mp], xp[mp])[0, 1] for i in range(len(Np))])
                    rec.update(null_mean=np.nanmean(nr), null_sd=np.nanstd(nr), null_p=(1 + (np.abs(nr) >= abs(r)).sum()) / (1 + len(nr)),
                               nullp_mean=np.nanmean(npr), nullp_sd=np.nanstd(npr))
                if f == 'prot':
                    # other kinases' protein
                    ro = np.array([pcor(ar, PKr[j])[0] for j in range(len(PKr))])
                    rpo = np.array([pcor(ap, PKr_p[j])[0] for j in range(len(PKr))])
                    own = (kin_m.uniprot.values == acc)
                    ro_o, rpo_o = ro[~own], rpo[~own]
                    rec.update(pct_other=np.nanmean(ro_o < r), pctp_other=np.nanmean(rpo_o < rp),
                               best_other=kin_m.kinase.values[~own][np.nanargmax(ro_o)], best_other_r=np.nanmax(ro_o),
                               best_otherp=kin_m.kinase.values[~own][np.nanargmax(rpo_o)], best_otherp_r=np.nanmax(rpo_o),
                               other_median_r=np.nanmedian(ro_o))
                rows.append(rec)
            # joint R2 (ranked): rna + prot (+ canonical site padj)
            cs = canon.get(k)
            X = {n_: feats.get(n_) for n_ in ['rna', 'prot']}
            if cs and f'{cs}|padj' in feats: X['site_padj'] = feats[f'{cs}|padj']
            X = {kk: rank(v) for kk, v in X.items() if v is not None}
            if X:
                M = np.column_stack(list(X.values())); m = np.isfinite(M).all(1) & np.isfinite(ar)
                if m.sum() >= 25:
                    Xd = np.column_stack([np.ones(m.sum()), M[m]]); y = ar[m]
                    b, *_ = np.linalg.lstsq(Xd, y, rcond=None); r2 = 1 - ((y - Xd @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
                    p_ = M.shape[1]; adj = 1 - (1 - r2) * (m.sum() - 1) / (m.sum() - p_ - 1)
                    jrows.append(dict(source=src, cohort=c, kinase=k, predictors='+'.join(X), n=int(m.sum()), R2=r2, R2_adj=adj))
        print(src, c, len(rows), flush=True)
pd.DataFrame(rows).to_parquet(f'{OUT}/cache/coupling_long.parquet')
pd.DataFrame(jrows).to_parquet(f'{OUT}/cache/joint_r2.parquet')
