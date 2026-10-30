#!/usr/bin/env python3
"""Per-sample kinase activity for CPTAC cohorts.

Input statistics (site x sample), built per cohort from phospho.parquet (log2), sites observed in >= 25% of the
cohort's tumours (and >= 10 tumours):
  tc         : x - median_tumours(x)  per site (tumours AND normals, both relative to the cohort tumour median)
  paired     : x_T - x_N for cases with tumour + normal adjacent phospho (no centring; raw paired difference)
  tc_padj    : residual of phospho ~ a + b*protein (OLS per site, fitted on tumours, protein = same UniProt
               accession, >= 10 complete tumour pairs); applied to tumours and normals with the tumour fit.
               Sites without protein measurement are dropped.
  paired_padj: tc_padj(T) - tc_padj(N) = dPhospho - b*dProtein for paired cases.
Methods:
  corev3 : kinaseActivityMain mvp core-v3 (edges = PSP + OmniPath curated + r2np predicted) -> z, z_core, n_substrates
  meanz  : per sample robust z of sites ((x - median)/1.4826 MAD over all sites in the sample); kinase score = mean z of
           its PSP-curated substrate sites (edges_psp_only, source 'cur'); NaN if < 5 measured substrates -> meanz, n
Outputs: /workspace/results/infra_cptac/activity/<variant>/<method>_{z,n,...}.parquet  (kinase x sample, all cohorts)
         and inputs  /workspace/results/infra_cptac/cohorts/<C>/input_<variant>.parquet (site x sample, float32)
"""
import os, sys, json, time
import numpy as np, pandas as pd

os.environ.setdefault('OMP_NUM_THREADS', '8')
MVP = '<HACKATHON_ROOT>/kinaseActivityMain/mvp'
sys.path.insert(0, MVP)
from kinact_mvp import KinaseActivityModel  # noqa
import torch
torch.set_num_threads(8)

OUT = '/workspace/results/infra_cptac'
COH = ['BRCA', 'CCRCC', 'COAD', 'GBM', 'HNSCC', 'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC']
MIN_FRAC, MIN_N, MIN_PAIRS, MIN_SUB = 0.25, 10, 10, 5


def build_inputs(c):
    d = f'{OUT}/cohorts/{c}'
    ph = pd.read_parquet(f'{d}/phospho.parquet').astype(np.float64)
    pr = pd.read_parquet(f'{d}/protein.parquet').astype(np.float64)
    T = [x for x in ph.columns if x.endswith('_T')]; N = [x for x in ph.columns if x.endswith('_N')]
    nT = ph[T].notna().sum(axis=1)
    ph = ph[(nT >= MIN_FRAC * len(T)) & (nT >= MIN_N)]
    med = ph[T].median(axis=1)
    tc = ph.sub(med, axis=0)
    out = {'tc': tc}
    paired_cases = sorted(set(x[:-2] for x in T) & set(x[:-2] for x in N))
    if paired_cases:
        out['paired'] = pd.DataFrame(ph[[p + '_T' for p in paired_cases]].values - ph[[p + '_N' for p in paired_cases]].values,
                                     index=ph.index, columns=[p + '_T' for p in paired_cases])
    # protein adjustment
    acc = ph.index.str.split('_').str[0]
    keep = acc.isin(pr.index)
    Y = ph[keep]; X = pr.reindex(acc[keep]).reindex(columns=Y.columns); X.index = Y.index
    yT, xT = Y[T].values, X[T].values
    m = np.isfinite(yT) & np.isfinite(xT)
    n = m.sum(1)
    with np.errstate(all='ignore'):
        my = np.where(m, yT, 0).sum(1) / n; mx = np.where(m, xT, 0).sum(1) / n
        dx = np.where(m, xT - mx[:, None], 0); dy = np.where(m, yT - my[:, None], 0)
        b = (dx * dy).sum(1) / (dx * dx).sum(1)
        a = my - b * mx
    ok = (n >= MIN_PAIRS) & np.isfinite(b)
    R = Y.values - (a[:, None] + b[:, None] * X.values)
    padj = pd.DataFrame(R[ok], index=Y.index[ok], columns=Y.columns)
    out['tc_padj'] = padj
    if paired_cases:
        out['paired_padj'] = pd.DataFrame(padj[[p + '_T' for p in paired_cases]].values - padj[[p + '_N' for p in paired_cases]].values,
                                          index=padj.index, columns=[p + '_T' for p in paired_cases])
    info = {'sites_filtered': int(len(ph)), 'sites_padj': int(ok.sum()), 'paired_cases': len(paired_cases),
            'median_slope_b': float(np.nanmedian(b[ok])), 'n_tumour': len(T), 'n_normal': len(N)}
    for k, v in out.items():
        v.astype(np.float32).to_parquet(f'{d}/input_{k}.parquet')
    return out, info


def meanz(V, ks):
    """V: site x sample DataFrame; ks: curated PSP edges (site, kinase). Returns (mean z, n) kinase x sample."""
    A = V.values
    med = np.nanmedian(A, axis=0); mad = 1.4826 * np.nanmedian(np.abs(A - med), axis=0)
    Z = (A - med) / mad
    e = ks[ks.site.isin(V.index)]
    idx = pd.Index(V.index)
    kin = sorted(e.kinase.unique())
    zs = np.full((len(kin), A.shape[1]), np.nan); ns = np.zeros((len(kin), A.shape[1]), int)
    for i, (k, g) in enumerate(e.groupby('kinase')):
        rows = idx.get_indexer(g.site.unique())
        sub = Z[rows]
        cnt = np.isfinite(sub).sum(0)
        with np.errstate(all='ignore'):
            zs[i] = np.where(cnt >= MIN_SUB, np.nansum(sub, 0) / cnt, np.nan)
        ns[i] = cnt
    return (pd.DataFrame(zs, index=kin, columns=V.columns), pd.DataFrame(ns, index=kin, columns=V.columns))


def main(cohorts):
    model = KinaseActivityModel.from_config(f'{OUT}/model_corev3/config.json')
    ks = pd.read_csv(f'{OUT}/edges/edges_psp_only.tsv.gz', sep='\t')
    ks = ks[ks.source.eq('cur')][['site', 'kinase']].drop_duplicates()
    res = {}; infos = {}
    for c in cohorts:
        t0 = time.time()
        ins, info = build_inputs(c)
        for var, V in ins.items():
            r = model.score(V, n_jobs=8)
            info[f'{var}_matched_lib_sites'] = r.info['n_library_sites_matched']
            mz, mn = meanz(V, ks)
            for name, df in (('corev3_z', r.z.T), ('corev3_zcore', r.z_core.T), ('corev3_n', r.n_substrates.T),
                             ('meanz_z', mz), ('meanz_n', mn)):
                res.setdefault((var, name), []).append(df)
        info['seconds'] = round(time.time() - t0, 1)
        infos[c] = info
        print(c, info, flush=True)
    for (var, name), dfs in res.items():
        od = f'{OUT}/activity/{var}'; os.makedirs(od, exist_ok=True)
        df = pd.concat(dfs, axis=1)
        df.index.name = 'kinase'
        (df.fillna(0).astype(np.int32) if name.endswith('_n') else df.astype(np.float32)).to_parquet(f'{od}/{name}.parquet')
    json.dump(infos, open(f'{OUT}/activity/activity_info.json', 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1:] or COH)
