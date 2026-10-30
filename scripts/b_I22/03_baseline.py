"""Baseline profiles of compendium cell lines: DepMap 24Q2 RNA (OmicsExpressionProteinCodingGenesTPMLogp1), CCLE proteomics (Nusinow 2020),
HPA cell-line RNA (fallback, HEK293/RPE1), Chronos 26Q1. Restricted to the compendium cell lines; all genes kept (later subset).
Run with python -I; inputs are untrusted downloads."""
import sys, zipfile, numpy as np, pandas as pd
DL = sys.argv[1]; OUT = sys.argv[2]
P = pd.read_csv(OUT + 'profiles_meta.tsv', sep='\t'); cells = sorted(P.cell.unique())
M = pd.read_csv(DL + 'Model.csv')
alias = {'P31/FUJ': 'P31FUJ', 'NTERA2': 'NTERA2CLD1', 'WIDR': 'HT29'}
s2id = dict(zip(M.StrippedCellLineName, M.ModelID))
cmap = {c: s2id.get(alias.get(c, c)) for c in cells}; cmap = {k: v for k, v in cmap.items() if v}
pd.Series(cmap).rename('ModelID').to_csv(OUT + 'cell_map.tsv', sep='\t')
# RNA
ex = pd.read_csv(DL + 'expr.csv', index_col=0)
ex = ex[ex.index.isin(cmap.values())]
ex.columns = [c.split(' (')[0] for c in ex.columns]
ex = ex.loc[:, ~ex.columns.duplicated()]
id2c = {}
for c, i in cmap.items(): id2c.setdefault(i, []).append(c)
rows = {c: ex.loc[i] for i in ex.index for c in id2c[i]}
R = pd.DataFrame(rows).T; R.to_parquet(OUT + 'rna_depmap.parquet'); print('RNA', R.shape)
# proteomics (Nusinow)
pr = pd.read_csv(DL + 'prot.csv.gz')
cc = [c for c in pr.columns if '_TenPx' in c and not c.endswith('_Peptides')]
ccle = dict(zip(M.ModelID, M.CCLEName))
want = {}
for c, i in cmap.items():
    nm = ccle.get(i)
    if isinstance(nm, str):
        hit = [x for x in cc if x.rsplit('_TenPx', 1)[0] == nm]
        if hit: want[c] = hit
G = pr.Gene_Symbol.values
prot = pd.DataFrame({c: pr[h].mean(1).values for c, h in want.items()}, index=G)
prot = prot[~pd.isna(prot.index)].groupby(level=0).mean().T
prot.to_parquet(OUT + 'prot_ccle.parquet'); print('PROT', prot.shape, sorted(want))
# HPA
z = zipfile.ZipFile(sys.argv[3]); h = pd.read_csv(z.open(z.namelist()[0]), sep='\t')
print(h.columns.tolist(), h.iloc[:2].to_string())
cl = h['Cell line'].unique()
hmap = {'HEK293': 'HEK293', 'RPE1': 'hTERT-RPE1', 'HELA': 'HeLa', 'MCF7': 'MCF7', 'HL60': 'HL-60', 'NTERA2': 'NTERA-2', 'A549': 'A549', 'U2OS': 'U-2 OS', 'A431': 'A-431', 'HCT116': 'HCT 116', 'K562': 'K-562', 'JURKAT': 'Jurkat'}
print([x for x in cl if any(k.lower() in x.lower() for k in ['293', 'rpe', 'hela', 'ntera', 'hl-60'])])
val = 'nTPM' if 'nTPM' in h.columns else h.columns[-1]
H = h[h['Cell line'].isin(hmap.values())].pivot_table(index='Cell line', columns='Gene name', values=val)
inv = {v: k for k, v in hmap.items()}; H.index = [inv[i] for i in H.index]
np.log2(H + 1).to_parquet(OUT + 'rna_hpa.parquet'); print('HPA', H.shape)
# Chronos
ch = pd.read_csv(sys.argv[4], index_col=0)
ch = ch[ch.index.isin(cmap.values())]; ch.columns = [c.split(' (')[0] for c in ch.columns]
C = pd.DataFrame({c: ch.loc[i] for i in ch.index for c in id2c[i]}).T; C.to_parquet(OUT + 'chronos.parquet'); print('CHRONOS', C.shape)
