"""Kinase Library (Johnson 2023 / Yaron-Barir 2024) PSSM log2 scores for every CPTAC site window.
Output: cache/motif_pct.parquet  (site_key x kinase gene; percentile within CPTAC sites of same residue class)
Score = sum_pos log2(PSSM) (+ log2 S/T favorability for ST kinases); lowercase/unknown ignored.
"""
import numpy as np, pandas as pd
D = '<HACKATHON_ROOT>/kinaseSubMain/external/kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases'
OUT = '/workspace/results/b_I20/cache'
info = pd.read_csv(f'{D}/kinase_data/kinome_information.tsv', sep='\t')
sm = pd.read_parquet('/workspace/results/infra_cptac/mapping/site_map.parquet')
sm = sm.dropna(subset=['site_key']).drop_duplicates('site_key')
win = sm.set_index('site_key').window15.str.upper()
res = win.str[7]
out = []
for lib, resset in [('ser_thr', 'ST'), ('tyrosine', 'Y')]:
    M = pd.read_csv(f'{D}/matrices/{lib}/{lib}_all_log2_matrices.tsv', sep='\t', index_col=0)
    w = win[res.isin(list(resset))]
    pos = sorted({int(c[:-1]) for c in M.columns})
    S = np.zeros((len(w), len(M)))
    for p in pos:
        aa = w.str[7 + p]
        for a in set(aa):
            col = f'{p}{a}'
            if col in M.columns:
                idx = (aa == a).values
                S[idx] += M[col].values[None, :]
    if lib == 'ser_thr':
        fav = pd.read_csv(f'{D}/matrices/ser_thr/st_favorability.tsv', sep='\t', index_col=0).loc[M.index]
        r = w.str[7].values
        S[r == 'S'] += np.log2(fav['s'].values.clip(1e-3))[None, :]
        S[r == 'T'] += np.log2(fav['t'].values.clip(1e-3))[None, :]
    df = pd.DataFrame(S, index=w.index, columns=M.index)
    pct = df.rank(pct=True) * 100
    out.append(pct.astype(np.float32))
pct = pd.concat(out, axis=0)  # ST kinases NaN on Y sites and vice versa
nm = info.set_index('MATRIX_NAME').GENE_NAME.to_dict()
pct.columns = [nm.get(c, c) for c in pct.columns]
pct = pct.T.groupby(level=0).max().T
pct.to_parquet(f'{OUT}/motif_pct.parquet')
print(pct.shape, pct.columns[:10].tolist())
