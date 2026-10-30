"""Task 1: instrument strength. protein ~ CNV(+purity+CIN) per gene per cohort; F = t^2.
Also: RNA ~ CNV, protein ~ focal CNV (gene minus arm median), protein ~ truncation (pooled).
Pooled within-cohort = fixed-effect IVW meta of per-cohort betas."""
import numpy as np, pandas as pd
from lib import *

kl = pd.read_csv('<HACKATHON_ROOT>/kinaseSubMain/results/labels/kinase_list.tsv', sep='\t')
kin = set(kl[kl.include.astype(str) == 'True'].uniprot_acc)
ann = pd.read_parquet(f'{CACHE}/gene_annot.parquet').drop_duplicates('uniprot').set_index('uniprot')
rows = []
for C in COHORTS:
    prot, cnv, foc, rna, tr, cov = [load(C, x) for x in ['prot', 'cnv', 'cnv_focal', 'rna', 'trunc', 'cov']]
    cv = cov.values
    G = prot.index.intersection(cnv.index)
    G = G[prot.loc[G].notna().sum(1) >= 30]
    Y = prot.loc[G].values.astype(float)
    out = pd.DataFrame(index=G)
    for nm, Z in [('cnv', cnv), ('focal', foc)]:
        b, t, n = rowwise_ols_t(Y, Z.loc[G].values.astype(float), cv)
        out[f'b_{nm}'] = b; out[f't_{nm}'] = t; out['n'] = n
        out[f'sd_{nm}'] = Z.loc[G].std(1).values
    Gr = G.intersection(rna.index)
    b, t, n = rowwise_ols_t(rna.loc[Gr].values.astype(float), cnv.loc[Gr].values.astype(float), cv)
    out.loc[Gr, 'b_rna_cnv'] = b; out.loc[Gr, 't_rna_cnv'] = t
    b, t, n = rowwise_ols_t(Y[np.isin(G, Gr)], rna.loc[G.intersection(Gr)].reindex(G[np.isin(G, Gr)]).values.astype(float), cv)
    out.loc[G[np.isin(G, Gr)], 't_prot_rna'] = t
    Gt = G.intersection(tr.index[tr.sum(1) >= 2])
    if len(Gt):
        b, t, n = rowwise_ols_t(prot.loc[Gt].values.astype(float), tr.loc[Gt].values.astype(float), cv)
        out.loc[Gt, 'b_trunc'] = b; out.loc[Gt, 't_trunc'] = t
        out.loc[Gt, 'n_trunc'] = tr.loc[Gt].sum(1).values
    out['cohort'] = C
    rows.append(out.rename_axis('uniprot').reset_index())
fs = pd.concat(rows)
fs['F_cnv'] = fs.t_cnv ** 2; fs['F_focal'] = fs.t_focal ** 2
fs['is_kinase'] = fs.uniprot.isin(kin)
fs['gene'] = fs.uniprot.map(ann.symbol); fs['arm'] = fs.uniprot.map(ann.arm)
fs.to_parquet(f'{RES}/tables/T1_first_stage_all_genes.parquet')


def meta(df, b, t):
    se = df[b] / df[t]
    w = 1 / se ** 2
    bm = (df[b] * w).sum() / w.sum()
    return pd.Series({'beta': bm, 'z': bm * np.sqrt(w.sum()), 'k': df[b].notna().sum()})


ks = fs[fs.is_kinase]
pooled = ks.groupby('uniprot').apply(lambda d: meta(d.dropna(subset=['t_cnv']), 'b_cnv', 't_cnv'))
pooled['F_pooled'] = pooled.z ** 2
tp = ks.dropna(subset=['t_trunc'])
tp = tp.groupby('uniprot').apply(lambda d: meta(d, 'b_trunc', 't_trunc')).rename(columns={'beta': 'beta_trunc', 'z': 'z_trunc', 'k': 'k_trunc'})
summ = ks.groupby('uniprot').agg(gene=('gene', 'first'), arm=('arm', 'first'), n_cohorts=('cohort', 'size'),
                                 n_strong=('F_cnv', lambda x: (x > 10).sum()), n_strong_focal=('F_focal', lambda x: (x > 10).sum()),
                                 med_F=('F_cnv', 'median'), med_b=('b_cnv', 'median'), med_t_rna=('t_rna_cnv', 'median'),
                                 med_t_prot_rna=('t_prot_rna', 'median'), n_trunc=('n_trunc', 'sum'))
summ = summ.join(pooled).join(tp).sort_values('F_pooled', ascending=False)
summ.to_csv(f'{RES}/tables/T1_first_stage_kinases.tsv', sep='\t')
print('kinases with protein+CNV:', len(summ))
print('kinase-cohort pairs F>10:', (ks.F_cnv > 10).sum(), 'of', len(ks), '; kinases strong in >=1 cohort:', (summ.n_strong >= 1).sum(),
      '>=2:', (summ.n_strong >= 2).sum(), '; pooled F>10:', (summ.F_pooled > 10).sum())
print('all genes F>10 frac:', (fs.F_cnv > 10).mean(), ' kinases:', (ks.F_cnv > 10).mean())
print('focal F>10 kinase-cohort pairs:', (ks.F_focal > 10).sum())
print('trunc: kinases with z_trunc<-3.16:', (summ.z_trunc < -np.sqrt(10)).sum(), summ[summ.z_trunc < -np.sqrt(10)].gene.tolist())
print('RNA check: median t rna~cnv among strong pairs', ks[ks.F_cnv > 10].t_rna_cnv.median(), ' weak', ks[ks.F_cnv <= 10].t_rna_cnv.median())
print(summ.head(40).to_string())
