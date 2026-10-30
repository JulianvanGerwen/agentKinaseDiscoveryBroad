#!/usr/bin/env python3
"""I08 step 1: kinase activity recomputed with ALL sites on protein kinases masked (so a kinase's activity never
uses its own (or a paralog's) phosphosites -> no circularity with own A-loop / own protein).
Methods: meanz (mean robust z of PSP-curated substrates, >=5) and corev3 (same model as infra_cptac).
Inputs: infra_cptac input_tc and input_tc_padj, tumours only.
Also: permutation null for meanz -- 200 random site sets per kinase, same size and same pY/pST composition.
Outputs: results/b_I08/cache/act_<method>_<var>.parquet (kinase x sample), null_<var>_<cohort>.npz
"""
import os, sys, numpy as np, pandas as pd
os.environ['OMP_NUM_THREADS'] = '8'; os.environ['MKL_NUM_THREADS'] = '8'
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
MVP = '<HACKATHON_ROOT>/kinaseActivityMain/mvp'; sys.path.insert(0, MVP)
OUT = '/workspace/results/b_I08/cache'; R = cl.ROOT
NPERM = 200; MIN_SUB = 5
kin = pd.read_csv(f'{R}/model_corev3/kinases.tsv', sep='\t')
kin_acc = set(kin.uniprot)
ks = pd.read_csv(f'{R}/edges/edges_psp_only.tsv.gz', sep='\t'); ks = ks[ks.source.eq('cur')][['site', 'kinase']].drop_duplicates()


def robz(V):
    A = V.values.astype(np.float64)
    med = np.nanmedian(A, 0); mad = 1.4826 * np.nanmedian(np.abs(A - med), 0)
    return (A - med) / mad


def setmean(Z, rows):
    sub = Z[rows]; cnt = np.isfinite(sub).sum(0)
    with np.errstate(all='ignore'):
        return np.where(cnt >= MIN_SUB, np.nansum(sub, 0) / cnt, np.nan)


def main():
    import torch; torch.set_num_threads(8)
    from kinact_mvp import KinaseActivityModel
    model = KinaseActivityModel.from_config(f'{R}/model_corev3/config.json')
    rng = np.random.default_rng(0)
    res = {}
    for var in ['tc_padj', 'tc']:
        for c in cl.COHORTS:
            V = cl.matrix(c, f'input_{var}')
            T = [x for x in V.columns if x.endswith('_T')]
            V = V[T]
            acc = V.index.str.split('_').str[0]
            V = V[~acc.isin(kin_acc)]          # mask all kinase-protein sites
            r = model.score(V, n_jobs=8)
            res.setdefault(('corev3', var), []).append(r.z.T)
            Z = robz(V); idx = pd.Index(V.index); isY = np.asarray(V.index.str.contains('_Y'))
            poolY = np.where(isY)[0]; poolST = np.where(~isY)[0]
            e = ks[ks.site.isin(V.index)]
            names, acts, nulls = [], [], []
            for k, g in e.groupby('kinase'):
                rows = idx.get_indexer(g.site.unique())
                a = setmean(Z, rows)
                if np.isfinite(a).sum() < 20: continue
                nY = isY[rows].sum(); nST = len(rows) - nY
                nl = np.empty((NPERM, Z.shape[1]), np.float32)
                for p in range(NPERM):
                    rr = np.concatenate([rng.choice(poolY, nY, replace=False), rng.choice(poolST, nST, replace=False)])
                    nl[p] = setmean(Z, rr)
                names.append(k); acts.append(a); nulls.append(nl)
            res.setdefault(('meanz', var), []).append(pd.DataFrame(np.array(acts), index=names, columns=V.columns))
            np.savez_compressed(f'{OUT}/null_{var}_{c}.npz', kin=np.array(names), null=np.array(nulls), samples=np.array(V.columns))
            print(var, c, V.shape, len(names), flush=True)
    for (m, var), dfs in res.items():
        df = pd.concat(dfs, axis=1); df.index.name = 'kinase'
        df.astype(np.float32).to_parquet(f'{OUT}/act_{m}_{var}.parquet')


if __name__ == '__main__':
    main()
