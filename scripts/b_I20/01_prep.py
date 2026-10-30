"""Per-cohort aligned tumour matrices for the MR/IV analysis.
cache/<C>.npz-like parquet set:
  prot.parquet   UniProt x case (tumours with phospho+protein+CNV)
  padj.parquet   site x case: phospho residualised on own protein + purity + CIN (OLS per site, n>=30)
  cnv.parquet    UniProt x case (gene CNV log2, mapped symbol->UniProt), cnv_focal (gene - arm median)
  rna.parquet    UniProt x case
  trunc.parquet  UniProt x case (0/1 truncating mutation)
  cov.parquet    case x [purity, cin]
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
OUT = '/workspace/results/b_I20/cache'
H = pd.read_csv('<HACKATHON_ROOT>/kinaseSubMain/external/hgnc/hgnc_complete_set.txt',
                sep='\t', usecols=['symbol', 'uniprot_ids', 'location', 'locus_group'], low_memory=False)
H = H.dropna(subset=['uniprot_ids'])
H['uniprot'] = H.uniprot_ids.str.split('|').str[0]
H['arm'] = H.location.str.extract(r'^(\d+|X|Y)([pq])').apply(lambda r: f'{r[0]}{r[1]}' if pd.notna(r[0]) else np.nan, axis=1)
H['chrom'] = H.location.str.extract(r'^(\d+|X|Y)')[0]
H[['symbol', 'uniprot', 'location', 'arm', 'chrom', 'locus_group']].to_parquet(f'{OUT}/gene_annot.parquet')
sym2u = H.drop_duplicates('symbol').set_index('symbol').uniprot
TRUNC = {'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site'}


def to_uni(df):
    df = df.copy(); df.index = df.index.map(sym2u)
    df = df[df.index.notna()]
    return df.groupby(level=0).mean()


for C in cl.COHORTS:
    ph = cl.matrix(C, 'phospho'); pr = cl.matrix(C, 'protein')
    cn = cl.matrix(C, 'cnv_log2'); rn = cl.matrix(C, 'rna')
    tum = [s for s in cl.tumour_ids(C) if s in pr.columns]
    case = {s: s[:-2] for s in tum}
    tum = [s for s in tum if case[s] in cn.columns]
    cases = [case[s] for s in tum]
    ph = ph[tum]; ph.columns = cases
    pr = pr[tum]; pr.columns = cases
    cn = cn[cases]
    rn = rn[[c for c in rn.columns if c in cases or c in tum]]
    rn.columns = [c[:-2] if c.endswith('_T') else c for c in rn.columns]
    rn = rn.reindex(columns=cases)
    clin = cl.table(C, 'clinical').drop_duplicates('case_id').set_index('case_id').reindex(cases)
    pur = clin.WES_purity.astype(float)
    pur = pur.fillna(pur.median())
    cin = cn.abs().mean(0)  # genome-wide CNV burden
    cov = pd.DataFrame({'purity': pur.values, 'cin': cin.values}, index=cases)
    # arm-level focal component
    g2arm = H.drop_duplicates('symbol').set_index('symbol').arm
    arm = cn.index.map(g2arm)
    armmed = cn.groupby(arm).median()
    foc = cn - armmed.reindex(arm).values
    foc.index = cn.index
    cnu, focu, rnu = to_uni(cn), to_uni(foc), to_uni(rn)
    # truncations
    mu = cl.table(C, 'mutations')
    mu = mu[mu.Variant_Classification.isin(TRUNC) & mu.case_id.isin(cases)]
    tr = pd.crosstab(mu.gene, mu.case_id).clip(upper=1).reindex(columns=cases, fill_value=0)
    tru = to_uni(tr).clip(upper=1)
    # protein-adjusted phospho, per site OLS on [1, own protein, purity, cin]
    uni = ph.index.str.split('_').str[0]
    keep = ph.notna().sum(1) >= 30
    keep &= uni.isin(pr.index)
    ph = ph[keep]; uni = uni[keep.values]
    P = pr.reindex(uni).values
    Y = ph.values.astype(np.float64)
    M = np.isfinite(Y) & np.isfinite(P)
    n = len(cases)
    X = np.stack([np.ones((len(Y), n)), np.nan_to_num(P), np.broadcast_to(cov.purity.values, (len(Y), n)),
                  np.broadcast_to(cov.cin.values, (len(Y), n))], axis=2)  # S x n x 4
    Xm = X * M[:, :, None]; Y0 = np.where(M, Y, 0)
    XtX = np.einsum('snj,snk->sjk', Xm, X) + 1e-8 * np.eye(4)
    Xty = np.einsum('snj,sn->sj', Xm, Y0)
    B = np.linalg.solve(XtX, Xty[:, :, None])[:, :, 0]
    R = Y - np.einsum('snj,sj->sn', X, B)
    R[~M] = np.nan
    okn = M.sum(1) >= 30
    padj = pd.DataFrame(R[okn].astype(np.float32), index=ph.index[okn], columns=cases)
    os.makedirs(f'{OUT}/{C}', exist_ok=True)
    padj.to_parquet(f'{OUT}/{C}/padj.parquet')
    pr.to_parquet(f'{OUT}/{C}/prot.parquet')
    cnu.to_parquet(f'{OUT}/{C}/cnv.parquet'); focu.to_parquet(f'{OUT}/{C}/cnv_focal.parquet')
    rnu.to_parquet(f'{OUT}/{C}/rna.parquet'); tru.to_parquet(f'{OUT}/{C}/trunc.parquet')
    cov.to_parquet(f'{OUT}/{C}/cov.parquet')
    print(C, 'n', n, 'sites', padj.shape[0], 'prot', pr.shape[0], 'trunc genes', tru.shape[0], flush=True)
