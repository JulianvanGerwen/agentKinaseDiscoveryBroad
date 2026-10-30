"""Cantley kinase-library (Johnson 2023 / Yaron-Barir 2024) motif scores for every compendium site.
score = sum_{pos -5..+4, !=0} log2 PSSM(pos, aa) (+ log2 S/T favourability for S/T kinases), computed from the
kinase_library 1.8.0 matrices (read-only). Percentile = rank among compendium sites of the same residue class
(S/T vs Y) -> stands in for the atlas percentile against the reference phosphoproteome.
Output: results/b_I05/cache/motif_pct.parquet (site x kinase GENE_NAME, uint8 percentile 0-100)."""
import numpy as np, pandas as pd
KL = '<HACKATHON_ROOT>/kinaseSubMain/external/kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases/'
OUT = '/workspace/results/b_I05/cache/'
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'residue', 'sequence_window'], low_memory=False)
S = S[S.sequence_window.str.len() == 15].reset_index(drop=True)
info = pd.read_csv(KL + 'kinase_data/kinome_information.tsv', sep='\t').set_index('MATRIX_NAME')
POS = [-5, -4, -3, -2, -1, 1, 2, 3, 4]
W = np.array([list(w.upper()) for w in S.sequence_window])  # n x 15, centre at 7
out = {}
for typ, res in [('ser_thr', ['S', 'T']), ('tyrosine', ['Y'])]:
    M = pd.read_csv(KL + f'matrices/{typ}/{typ}_all_log2_matrices.tsv', sep='\t', index_col=0)
    m = S.residue.isin(res).values
    Wm = W[m]; sc = np.zeros((m.sum(), len(M)), np.float32)
    for p in POS:
        aa = Wm[:, 7 + p]
        for a in np.unique(aa):
            col = f'{p}{a}'
            if col in M.columns:
                sc[aa == a] += M[col].values.astype(np.float32)
    if typ == 'ser_thr':
        fav = pd.read_csv(KL + 'matrices/ser_thr/st_favorability.tsv', sep='\t', index_col=0).reindex(M.index)
        r = np.asarray(S.residue.values[m], dtype=object)
        sc += np.where(r[:, None] == 'S', np.log2(fav['s'].values)[None, :], np.log2(fav['t'].values)[None, :]).astype(np.float32)
    pct = (pd.DataFrame(sc).rank(pct=True).values * 100).round().astype(np.uint8)
    genes = [info.loc[k, 'GENE_NAME'] if k in info.index else k for k in M.index]
    df = pd.DataFrame(pct, index=S.site.values[m], columns=genes)
    df = df.T.groupby(level=0).max().T  # duplicate gene names (rare)
    out[typ] = df
    print(typ, df.shape)
P = pd.concat(out.values(), axis=0)  # S/T kinases NaN for Y sites and vice versa
P = P.fillna(0).astype(np.uint8)
P.to_parquet(OUT + 'motif_pct.parquet')
print(P.shape, P.columns[:20].tolist())
