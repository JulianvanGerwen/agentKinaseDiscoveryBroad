"""I29 replication: Chen 2020 Cell TW LUAD (Table S1, log2 T/N). Reads small extracted tsvs only."""
import numpy as np, pandas as pd, statsmodels.api as sm, re
from scipy import stats
D = '/workspace/cache/b_I29/chen/'; OUT = '/workspace/results/b_I29/'
def z(s): s = pd.to_numeric(s, errors='coerce'); return (s - s.mean()) / s.std()
clin = pd.read_csv(D + 'TableS1A.tsv', sep='\t').set_index('ID')
snv = pd.read_csv(D + 'TableS1C.tsv', sep='\t', skiprows=1).set_index('Gene')
snv = snv[[c for c in snv.columns if re.match(r'P\d+', str(c))]]
rasmut = snv.reindex(['KRAS', 'NRAS', 'BRAF', 'HRAS', 'NF1']).notna().any(axis=0)
rna = pd.read_csv(D + 'TableS1D.tsv', sep='\t').set_index('gene'); rna = rna[[c for c in rna.columns if re.match(r'P\d+$', str(c))]]
prot = pd.read_csv(D + 'TableS1E.tsv', sep='\t').set_index('Gene'); prot = prot[[c for c in prot.columns if re.match(r'P\d+$', str(c))]]
G = pd.read_csv(D + 'TableS1G_acc.tsv', sep='\t').set_index('ID'); G = G[[c for c in G.columns if re.match(r'P\d+$', str(c))]].apply(pd.to_numeric, errors='coerce')
F = pd.read_csv(D + 'TableS1F_acc.tsv', sep='\t'); F.index = F['Modifications in Master Protein'] + ' | ' + F['Modifications'].str.contains('Y15').map({True: 'Y15', False: ''})
F = F[[c for c in F.columns if re.match(r'P\d+$', str(c))]].apply(pd.to_numeric, errors='coerce')
pts = sorted(set(clin.index))
A = pd.DataFrame(index=pts)
A['EGFR'] = clin.EGFR_Status.reindex(pts); A['hist'] = clin['Histology Type'].reindex(pts); A['smoke'] = clin['Smoking Status'].reindex(pts); A['stage'] = clin.Stage.reindex(pts)
A['ras'] = rasmut.reindex(pts).fillna(False).astype(bool)
A['group'] = np.where(A.EGFR != 'WT', 'EGFR', np.where(A.ras, 'RAS', 'WT'))
A['smoker'] = (A.smoke != 'Nonsmoke').astype(float)
A['stage_num'] = A.stage.astype(str).str.extract(r'(IV|III|II|I)')[0].map({'I': 1, 'II': 2, 'III': 3, 'IV': 4}).astype(float)
site = {'MAPK3_T202': 'P27361 1xPhospho [T202]', 'MAPK1_T185': 'P28482 1xPhospho [T185]', 'MAPK1_T190(=Y187?)': 'P28482 1xPhospho [T190]'}
for k, v in site.items(): A[k] = G.loc[v].reindex(pts) if v in G.index else np.nan
pep = {'MAPK3_pTpY_T202Y204': 'P27361 [T202; Y204] | Y15', 'MAPK1_pTpY_T185Y187': 'P28482 [T185; Y187] | Y15',
       'MAPK1_dbl_T185T190(ambig)': 'P28482 [T185; T190] | ', 'MAPK3_dbl_T202T207(ambig)': 'P27361 [T202; T207] | ', 'MAPK3_mono_T202_pep': 'P27361 [T202] | '}
for k, v in pep.items():
    m = F.loc[[i for i in F.index if i == v]]
    A[k] = m.mean(axis=0).reindex(pts) if len(m) else np.nan
for g in ['ETV4', 'ETV5', 'SPRED1', 'DUSP6', 'DUSP4', 'SPRY4', 'SPRY2', 'ETV1', 'EGR1', 'FOS', 'MAPK1', 'MAPK3', 'EGFR']:
    A['rna_' + g] = rna.loc[g].reindex(pts) if g in rna.index else np.nan
for g in ['MAPK1', 'MAPK3', 'DUSP6', 'EGFR', 'SPRED1']:
    A['prot_' + g] = prot.loc[g].reindex(pts) if g in prot.index else np.nan
A['ETVc'] = pd.concat([z(A['rna_' + g]) for g in ['ETV4', 'ETV5', 'SPRED1', 'DUSP6']], axis=1).mean(axis=1)
A['pERK'] = pd.concat([z(A.MAPK3_T202), z(A.MAPK1_T185)], axis=1).mean(axis=1)
# protein-adjusted
for s, p in [('MAPK3_T202', 'prot_MAPK3'), ('MAPK1_T185', 'prot_MAPK1')]:
    d = A[[s, p]].dropna(); b = np.polyfit(d[p], d[s], 1); A[s + '_padj'] = A[s] - (b[0] * A[p] + b[1])
A['pERK_padj'] = pd.concat([z(A.MAPK3_T202_padj), z(A.MAPK1_T185_padj)], axis=1).mean(axis=1)
A['pERK_dbl_ambig'] = pd.concat([z(A['MAPK1_dbl_T185T190(ambig)']), z(A['MAPK3_dbl_T202T207(ambig)'])], axis=1).mean(axis=1)
A.to_csv(OUT + 'tables/T5_chen_tumour_table.tsv', sep='\t')
print(A.groupby(['hist', 'group']).size())
print(A[[c for c in A.columns if 'MAPK' in c or c.startswith('pERK')]].notna().sum())
rows = []
def fit(y, mask, covs=(), groups=('EGFR', 'RAS'), label=''):
    d = A[mask & A.group.isin(list(groups) + ['WT'])][[y, 'group', *covs]].dropna()
    if (d.group == 'WT').sum() < 4 or len(d) < 12: return
    X = pd.DataFrame({g: (d.group == g) * 1.0 for g in groups if (d.group == g).sum() >= 3}, index=d.index)
    for c in covs: X[c] = z(d[c]).fillna(0)
    Y = (d[y] - A[y].mean()) / A[y].std()
    r = sm.OLS(Y, sm.add_constant(X)).fit(cov_type='HC3')
    for g in X.columns:
        if g in covs: continue
        rows.append(dict(outcome=y, model=label, contrast=f'{g}_vs_WT', beta=r.params[g], lo=r.conf_int().loc[g, 0], hi=r.conf_int().loc[g, 1], p=r.pvalues[g],
                         n_g=int((d.group == g).sum()), n_ref=int((d.group == 'WT').sum())))
OUTS = ['pERK', 'pERK_padj', 'MAPK3_T202', 'MAPK1_T185', 'MAPK1_T190(=Y187?)', 'MAPK3_T202_padj', 'MAPK1_T185_padj', 'MAPK3_pTpY_T202Y204', 'MAPK1_pTpY_T185Y187', 'pERK_dbl_ambig', 'MAPK3_mono_T202_pep',
        'prot_MAPK1', 'prot_MAPK3', 'prot_DUSP6', 'prot_EGFR', 'ETVc', 'rna_ETV4', 'rna_ETV5', 'rna_SPRED1', 'rna_DUSP6', 'rna_DUSP4', 'rna_SPRY4', 'rna_EGR1', 'rna_FOS', 'rna_EGFR']
adc = A['hist'] == 'ADC'
for o in OUTS:
    if A[o].notna().sum() < 12: continue
    fit(o, adc, label='ADC (primary)'); fit(o, adc, ['stage_num', 'smoker'], label='ADC +stage +smoking')
    fit(o, adc & (A.smoke == 'Nonsmoke'), label='ADC never-smokers')
    A['group2'] = A.group; A['group'] = np.where(A.group2 == 'EGFR', 'EGFR', 'WT'); fit(o, adc, groups=('EGFR',), label='ADC EGFR vs all non-EGFR (deviation: WT n=5)'); A['group'] = A.group2; fit(o, A['hist'].notna(), label='all histologies')
T = pd.DataFrame(rows)
for c in ['p']:
    T['p_one_sided_pred'] = np.where(T.outcome.str.contains('ETV|SPRED|DUSP6|rna_ETV'), np.where(T.beta > 0, T.p / 2, 1 - T.p / 2), np.where(T.beta < 0, T.p / 2, 1 - T.p / 2))
T.to_csv(OUT + 'tables/T6_replication_chen.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
print(T[T.contrast == 'EGFR_vs_WT'].round(3).to_string())
print(T[T.contrast == 'RAS_vs_WT'][['outcome', 'model', 'beta', 'p', 'n_g', 'n_ref']].round(3).to_string())
E = A[(A.group == 'EGFR') & adc]
for x in ['rna_DUSP6', 'prot_DUSP6']:
    for y in ['pERK', 'pERK_padj', 'rna_ETV4', 'ETVc']:
        d = E[[x, y]].dropna()
        if len(d) > 10: r, p = stats.spearmanr(d[x], d[y]); print('within EGFR ADC', x, y, round(r, 3), round(p, 3), len(d))
d = A[adc][['pERK', 'rna_ETV4']].dropna(); print('all ADC pERK vs ETV4', stats.spearmanr(d.pERK, d.rna_ETV4))
