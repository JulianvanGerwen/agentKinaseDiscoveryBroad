"""audit_I10 step 5: is the MAPK/EGF footprint axis a cold-ischaemia (pre-analytical) signature?
Ischaemia footprint = top-200 |mean robust z| sites over the 3 Mertins 2014 (PMID 24719451) cold-ischaemia experiments in the compendium
(sign-consistent in all measured, >=2 measured). Per cohort (tc_padj): own scores; r(MEK, ISCH), r(EGF, ISCH); site overlap and sign concordance;
MEK/EGF coherence after (a) dropping footprint sites in the ischaemia top-500 and (b) regressing the ischaemia score out of every site.
Output: A5_ischaemia.tsv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from common import *
W = pd.read_parquet('/workspace/cache/b_I10/fp_sites.parquet')
ex = [c for c in W.columns if c.startswith('24719451|') and 'ischaemia' in c]
I = W[ex]; n = I.notna().sum(1); mu = I.mean(1); cons = (np.sign(I).eq(np.sign(mu), axis=0) | I.isna()).all(1)
ok = (n >= 2) & cons
rank = mu[ok].abs().sort_values(ascending=False)
ISCH = pd.Series(np.sign(mu[rank.index[:200]]).values, index=rank.index[:200]); ISCH500 = set(rank.index[:500])
del W
rows = []
for g in ['KI_MEK', 'EGF', 'GF_all', 'MAPKaxis', 'TNF', 'DNA_damage', 'phosphatase_inh', 'metformin', 'GPCR_ligand']:
    s = fp_set(g); sh = s.index.intersection(ISCH.index)
    rows.append(dict(cohort='ALL', set=g, what='overlap with ischaemia top-200 (n, sign-concordance)', value=f'{len(sh)} ({(s[sh] == ISCH[sh]).mean() if len(sh) else np.nan:.2f})'))
    sh5 = s.index[s.index.isin(ISCH500)]; rows.append(dict(cohort='ALL', set=g, what='n footprint sites in ischaemia top-500', value=len(sh5)))
for coh in COH:
    Z = standardise(cohort(coh)); C0 = Coh(Z)
    si = score(Z, ISCH)
    # regress ischaemia score out of every site
    v = si.reindex(Z.columns); vv = (v - v.mean()) / v.std()
    X = Z.values; M = np.isfinite(X); Xr = np.where(M, X, np.nan)
    vz = np.nan_to_num(vv.values)
    with np.errstate(all='ignore'):
        beta = np.nansum(Xr * vz, 1) / (M * vz ** 2).sum(1)
    Zr = standardise(pd.DataFrame(Xr - beta[:, None] * vz, index=Z.index, columns=Z.columns)); Cr = Coh(Zr)
    for g in ['KI_MEK', 'EGF', 'GF_all', 'TNF', 'DNA_damage', 'phosphatase_inh', 'metformin', 'GPCR_ligand']:
        s = fp_set(g); s = s[s.index.isin(C0.idx)]
        sc = score(Z, s); r = pd.concat([sc, si], axis=1).dropna().corr().iloc[0, 1]
        ii = np.array([C0.idx[t] for t in s.index])
        s2 = s[~s.index.isin(ISCH500)]; ii2 = np.array([C0.idx[t] for t in s2.index])
        rows += [dict(cohort=coh, set=g, what='r(score, ischaemia score)', value=r),
                 dict(cohort=coh, set=g, what='coherence base', value=C0.coherence(ii, s.values)[0]),
                 dict(cohort=coh, set=g, what='coherence w/o ischaemia top-500 sites', value=C0.coherence(ii2, s2.values)[0]),
                 dict(cohort=coh, set=g, what='coherence after regressing ischaemia score', value=Cr.coherence(ii, s.values)[0])]
    sI = ISCH[ISCH.index.isin(C0.idx)]; rows.append(dict(cohort=coh, set='ISCHAEMIA', what='coherence base', value=C0.coherence(np.array([C0.idx[t] for t in sI.index]), sI.values)[0]))
    print(coh, flush=True)
T = pd.DataFrame(rows); T.to_csv(RA + 'A5_ischaemia.tsv', sep='\t', index=False)
t = T[T.cohort != 'ALL'].copy(); t['value'] = t.value.astype(float)
pd.set_option('display.width', 250); print(T[T.cohort == 'ALL'].to_string())
print(t.pivot_table(index='set', columns='what', values='value', aggfunc='median').round(3).to_string())
print(t[t.what.str.startswith('r(')].pivot_table(index='set', columns='cohort', values='value').round(2).to_string())
