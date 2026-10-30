"""Positive control for the H2 machinery: per unit (study x cell), fraction of a gene's compendium sites detected in that unit's experiments
vs DepMap RNA z of the gene (all genes with >=3 sites); matched vs permuted cell labels; also restricted to kinase genes only."""
import numpy as np, pandas as pd
from scipy.stats import rankdata
CA = '/workspace/cache/b_I22/'; OUT = '/workspace/results/b_I22/tables/'
Mk = pd.read_parquet(CA + 'measured_mask.parquet')
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'gene'], low_memory=False).set_index('site')
R = pd.read_parquet(CA + 'rna_depmap.parquet'); kst = list(pd.read_parquet(CA + 'profiles.parquet').columns)
meta = pd.DataFrame({'col': Mk.columns}); meta['kind'] = meta.col.str.split('::').str[0]; meta['exp'] = meta.col.str.split('::').str[1]
meta['study'] = meta.exp.str.split('|').str[0]; meta['cell'] = meta.exp.str.split('|').str[1]; meta['unit'] = meta.study + '|' + meta.cell
meta = meta[meta.cell.isin(R.index)]
g = S.gene.reindex(Mk.index).values
U = Mk[meta.col].T.groupby(meta.unit.values).max().T  # site measured in any exp of unit
D = U.groupby(g).mean(); D = D[(pd.Series(g).value_counts().reindex(D.index) >= 3).values]
um = meta.groupby('unit').cell.first()
cells = sorted(um.unique()); genes = [x for x in D.index if x in R.columns]
Rz = R.loc[cells, genes]; Rz = (Rz - Rz.mean()) / Rz.std()
Dv = D.loc[genes, um.index].T; Dv = Dv - Dv.mean()  # detection relative to other units
out = []
for lab, gs in [('all genes', genes), ('kinase genes (Cantley S/T)', [x for x in genes if x in kst])]:
    Dr = np.apply_along_axis(rankdata, 1, Dv[gs].values); Rr = {c: rankdata(Rz.loc[c, gs].fillna(0).values) for c in cells}
    stat = lambda a: np.nanmean([np.corrcoef(Dr[i], Rr[a[c]])[0, 1] for i, c in enumerate(um.values)])
    obs = stat({c: c for c in cells}); rng = np.random.default_rng(0)
    null = np.array([stat(dict(zip(cells, rng.permutation(cells)))) for _ in range(300)])
    out.append(dict(test=lab, n_genes=len(gs), n_units=len(um), mean_rho=obs, null_mean=null.mean(), null_sd=null.std(), perm_p=(np.sum(null >= obs) + 1) / 301))
T = pd.DataFrame(out); T.to_csv(OUT + 'T7_posctrl_detection_vs_rna.tsv', sep='\t', index=False); print(T.round(3).to_string())
