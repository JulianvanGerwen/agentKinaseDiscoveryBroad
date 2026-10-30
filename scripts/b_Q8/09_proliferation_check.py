"""b_Q8 step 9: is the tumour gain of DNA-damage / CDK2 / DDR coherence just proliferation? Within each group (paired, padj),
regress per-sample proliferation covariates (PSP CDK1 + CDK2 substrate scores and mitosis_cellcycle footprint score, each computed on
sites NOT in the tested set) out of every site, then recompute coherence T vs NAT with matched random-site nulls (bins coverage x SD,
100 draws) and 100 within-patient label swaps."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from common import *
rng = np.random.default_rng(3)
PS = psp_sets(15)
TEST = {'FP:DNA_damage': fp_set('DNA_damage'), 'FP:KI_DDR': fp_set('KI_DDR'), 'FP:replication_stress': fp_set('replication_stress'),
        'PSP:CDK2': pd.Series(1.0, index=PS['CDK2']), 'FP:KI_MEK': fp_set('KI_MEK'), 'FP:EGF': fp_set('EGF')}
PRO = {'CDK1': pd.Series(1.0, index=PS['CDK1']), 'CDK2': pd.Series(1.0, index=PS['CDK2']), 'mit': fp_set('mitosis_cellcycle')}
rows = []
for coh in MAINC:
    ph, pr, T, N, paired = load_raw(coh)
    cs = sorted(c for c in paired if f'{c}_N' in N and f'{c}_T' in T); n = len(cs); mp = max(10, n // 2)
    for nm, sg in TEST.items():
        out = {}
        for grp, smp in [('T', [c + '_T' for c in cs]), ('N', [c + '_N' for c in cs])]:
            Z0 = process(ph, pr, smp, padj=True, k=0, minobs=max(8, int(0.4 * n)))
            cov = pd.DataFrame({k: score(Z0, v[~v.index.isin(sg.index)]) for k, v in PRO.items() if not (nm == 'PSP:CDK2' and k == 'CDK2')}).fillna(0)
            Z = process(ph, pr, smp, padj=True, k=0, cov=cov, minobs=max(8, int(0.4 * n)))
            out[grp] = Z
        ix = out['T'].index.intersection(out['N'].index); ZT, ZN = out['T'].loc[ix], out['N'].loc[ix]
        CT, CN = Coh(ZT, mp), Coh(ZN, mp)
        s = sg[sg.index.isin(CT.idx)]; ii = np.array([CT.idx[x] for x in s.index]); sv = s.values
        cT, cN = CT.coherence(ii, sv)[0], CN.coherence(ii, sv)[0]
        q = lambda v, k: pd.qcut(pd.Series(v).rank(method='first'), k, labels=False).values
        b = q(ph.loc[ix, [c + '_T' for c in cs]].notna().mean(1).values, 4) * 4 + q(ph.loc[ix].std(1).values, 4)
        pool = {v: np.where(b == v)[0] for v in np.unique(b)}
        nT, nN = [], []
        for _ in range(100):
            jj = np.array([rng.choice(pool[v]) for v in b[ii]]); nT.append(CT.coherence(jj, sv)[0]); nN.append(CN.coherence(jj, sv)[0])
        perm = []
        for _ in range(100):
            sw = rng.random(n) < 0.5
            A = np.where(sw[None, :], ZN.values, ZT.values); B = np.where(sw[None, :], ZT.values, ZN.values)
            perm.append(Coh(pd.DataFrame(A, index=ix), mp).coherence(ii, sv)[0] - Coh(pd.DataFrame(B, index=ix), mp).coherence(ii, sv)[0])
        perm = np.array(perm)
        rows.append(dict(cohort=coh, set=nm, coh_T=cT, coh_N=cN, z_T=(cT - np.mean(nT)) / np.std(nT), z_N=(cN - np.mean(nN)) / np.std(nN),
                         delta=cT - cN, z_delta_perm=(cT - cN - perm.mean()) / perm.std()))
    print(coh, flush=True)
D = pd.DataFrame(rows); D.to_csv(TQ + 'T11_proliferation_adjusted_coherence.tsv', sep='\t', index=False)
g = D.groupby('set'); print(pd.DataFrame(dict(z_T=g.z_T.median(), z_N=g.z_N.median(), delta=g.delta.median(), meta_z=g.z_delta_perm.sum() / np.sqrt(g.size()),
                                             n_pos=g.z_delta_perm.apply(lambda x: (x > 1.96).sum()))).round(2))
