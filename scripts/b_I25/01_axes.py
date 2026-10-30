"""Reference axes (site scores + sets) and site strata for matching. Outputs cache/b_I25/axes.parquet."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = '/workspace/cache/b_I25/'
C = load.conditions()
AX = {
 'M': ['20068231_qPTM_2_log2', '18669648_Ochoa_267_72_log2', '18691976_qPTM_1_log2', '28674419_qPTM_2_log2', '29052541_qPTM_2_log2'],
 'MKI': ['26055452_qPTM_1_log2', '27134283_qPTM_2_log2', '26270265_qPTM_1_log2'] + [f'21712546_qPTM_{i}_log2' for i in (1, 2, 4, 5, 6)],
 'INTER': ['32051232_1_log2', '32051232_2_log2', '32051232_3_log2', '29052541_qPTM_1_log2', '29052541_qPTM_3_log2', '23322592_qPTM_1_log2', '18669648_Ochoa_72_72_log2', '29154753_qPTM_1_log2'],
 'MET24': C.conditionID[(C.study_pmid == '36781298') & (C.time_min == 1440)].tolist(),
 'MET30': C.conditionID[(C.study_pmid == '36781298') & (C.time_min == 30)].tolist(),
 'CULT3D': ['30061712_qPTM_1_log2', '30061712_qPTM_2_log2'],
 'SERUM': C.conditionID[C.study_pmid == '29970458'].iloc[:0].tolist() + ['29970458_1_log2', '29970458_2_log2', '29970458_3_log2', '29970458_4_log2'],
 'OIS': ['24961811_qPTM_2_log2', '24961811_qPTM_4_log2', '24961811_qPTM_5_log2', '24961811_qPTM_6_log2'],
}
SIGN = {'MKI': -1}
allc = sorted({c for v in AX.values() for c in v})
X = load.sites_matrix(allc)
Z = (X - X.median()) / (1.4826 * (X - X.median()).abs().median())
res = {}
for ax, cols in AX.items():
    z = Z[cols] * SIGN.get(ax, 1)
    n = z.notna().sum(1); m = z.mean(1)
    agree = np.where(m > 0, (z > 0).sum(1), (z < 0).sum(1)) / n.replace(0, np.nan)
    res[f'{ax}_score'] = m; res[f'{ax}_n'] = n; res[f'{ax}_agree'] = agree
A = pd.DataFrame(res)
def sett(ax, sgn=None, thr=1.5):
    m, ag = A[f'{ax}_score'], A[f'{ax}_agree']
    ok = (ag >= 0.75) | (A[f'{ax}_n'] == 1)
    if sgn == 1: return (m >= thr) & ok
    if sgn == -1: return (m <= -thr) & ok
    return (m.abs() >= thr) & ok
A['CC_Mup'] = sett('M', 1); A['CC_Mdown'] = sett('M', -1); A['CC_MKI'] = sett('MKI', 1); A['CC_inter'] = sett('INTER')
A['CC_union'] = A.CC_Mup | A.CC_Mdown | A.CC_MKI
for ax in ['MET24', 'MET30', 'CULT3D', 'SERUM', 'OIS']: A[f'S_{ax}'] = sett(ax)
A['S_MET24up'] = sett('MET24', 1); A['S_MET24down'] = sett('MET24', -1)
# universes
A['U_CC'] = (A.M_n > 0) | (A.MKI_n > 0)
for ax in ['M', 'MKI', 'INTER', 'MET24', 'MET30', 'CULT3D', 'SERUM', 'OIS']: A[f'U_{ax}'] = A[f'{ax}_n'] > 0
S = load.sites().set_index('site')
A = A.reindex(S.index)
for c in A.columns:
    if A[c].dtype == object or c.startswith(('CC_', 'S_', 'U_')): A[c] = A[c].fillna(False).astype(bool)
# motif energy set
MP = pd.read_parquet('/workspace/results/b_I05/cache/motif_pct.parquet')
ek = [k for k in ['PRKAA1', 'PRKAA2', 'MARK1', 'MARK2', 'MARK3', 'MARK4', 'SIK1', 'SIK2', 'SIK3', 'NUAK1', 'NUAK2', 'BRSK1', 'BRSK2'] if k in MP.columns]
print('energy motif kinases', ek)
A['S_AMPKmotif'] = (MP.reindex(A.index)[ek].max(1) >= 95).fillna(False).values
A['U_AMPKmotif'] = MP.reindex(A.index)[ek].notna().any(axis=1).values & S.residue.isin(['S', 'T']).values
HIPPO = 'YAP1 WWTR1 LATS1 LATS2 STK3 STK4 MOB1A MOB1B SAV1 NF2 AMOT AMOTL1 AMOTL2 WWC1 PTPN14 FRMD6 TEAD1 TEAD2 TEAD3 TEAD4'.split()
A['S_HIPPO'] = S.gene.isin(HIPPO).values; A['U_HIPPO'] = True
# strata: residue x coverage decile x variance decile (variance over non-test conditions)
M = pd.read_csv('/workspace/results/audit_I19/cache/inh_meta.tsv', sep='\t')
test = set(M.conditionID)
import pyarrow.parquet as pq
cols = [c for c in pq.read_schema('/workspace/results/infra_compendium/site_matrix.parquet').names if c != 'site' and c not in test]
sd = []; 
for i in range(0, len(cols), 300):
    Xi = load.sites_matrix(cols[i:i + 300]); sd.append(Xi)
Xn = pd.concat(sd, axis=1)
v = Xn.std(1, skipna=True); nn = Xn.notna().sum(1)
A['var_sd'] = v.reindex(A.index).values
A['ncond'] = S.n_conditions.values; A['residue'] = S.residue.values
A['cov_dec'] = pd.qcut(A.ncond.rank(method='first'), 10, labels=False)
vv = A.var_sd.fillna(A.var_sd.median())
A['var_dec'] = pd.qcut(vv.rank(method='first'), 10, labels=False)
A['stratum'] = A.residue.astype(str) + '_' + A.cov_dec.astype(str) + '_' + A.var_dec.astype(str)
A.index.name = 'site'
A.to_parquet(OUT + 'axes.parquet')
sets = [c for c in A.columns if c.startswith(('CC_', 'S_'))]
print(pd.DataFrame({'n_set': A[sets].sum(), }).T.to_string())
for c in ['U_CC', 'U_M', 'U_MKI', 'U_INTER', 'U_MET24', 'U_CULT3D', 'U_SERUM', 'U_OIS']: print(c, A[c].sum())
print(A.groupby('CC_Mup').var_sd.median())
# sanity: classical mitotic sites
for s in ['P06493_Y15', 'Q8TD19_S869', 'P16949_S16', 'P16949_S38', 'P53350_T210']:
    print(s, A.loc[s, ['M_score', 'M_n', 'MKI_score', 'CC_Mup', 'CC_Mdown', 'CC_MKI']].tolist())
