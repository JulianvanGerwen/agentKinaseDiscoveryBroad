"""audit_I10 step 1: independent recomputation of footprint / PSP-set coherence in CPTAC tumours + confound stress tests.
Data: input_tc_padj (phospho residualised on own protein, infra_cptac), tumours, sites >=25% coverage; site-wise z over observed values;
pairwise-complete Pearson (n_ij>=20). Signed coherence = mean over cross-protein pairs of s_i s_j r_ij.
Conditions: base | coab (drop pairs whose PROTEINS correlate |r|>0.3 in that cohort) | string (drop STRING physical>=700 pairs) |
            coab+string | pc1/pc3/pc5 (top-k sample PCs of all sites regressed out of every site) | pc5+coab+string.
Nulls (200 per set x cohort x condition), signs of the real set kept:
  nA = coverage x SD bins (as b_I10)    nB = coverage x SD x communality(top-5 PC R^2) x residue(Y vs S/T) x proline+1 bins
Also: random single-experiment 'footprints' (top-200 |z| sites of random compendium experiments) as a biology-agnostic reference.
Outputs: results/audit_I10/A1_coherence.tsv, A1_random_experiments.tsv, A1_site_features.tsv"""
import numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
from common import *
rng = np.random.default_rng(11); NN = 200
FPS = ['KI_MEK', 'EGF', 'GF_all', 'MAPKaxis', 'TNF', 'DNA_damage', 'phosphatase_inh', 'GPCR_ligand', 'metformin', 'KI_EGFR']
PSPN = ['CDK1', 'CDK2', 'MAPK1', 'MAPK3', 'AKT1', 'PRKACA', 'CSNK2A1', 'RPS6KA1', 'MAPKAPK2', 'SRC', 'MTOR', 'ATM', 'PRKCA', 'GSK3B', 'CHEK1']
SETS = {f'FP:{g}': fp_set(g) for g in FPS}
PSP = psp_sets(15)
for k, v in PSP.items(): SETS[f'PSP:{k}'] = v
NULL_SETS = set([f'FP:{g}' for g in FPS] + [f'PSP:{k}' for k in PSPN])
SA = site_annot()
# random compendium experiments
W = pd.read_parquet('/workspace/cache/b_I10/fp_sites.parquet')
rexp = rng.choice(W.columns, 150, replace=False)
REXP = {}
for e in rexp:
    z = W[e].dropna(); top = z.abs().sort_values(ascending=False).index[:200]
    REXP[e] = pd.Series(np.sign(z[top].values), index=top)
del W
STR = pd.read_parquet(CA + 'string_phys700.parquet')
rows, rrows, feats = [], [], []
for coh in COH:
    t0 = time.time()
    x = cohort(coh); Z = standardise(x)
    pr = cl.matrix(coh, 'protein'); pr = pr[[c for c in Z.columns if c in pr.columns]]
    Pc = Coh(standardise(pr[pr.notna().sum(1) >= 20]))
    # PCs and site features
    Rk = {}
    for k in [1, 3, 5]:
        Rk[k], F, s = residualise_pcs(Z, k)
        if k == 5: F5, s5 = F, s
    C0 = Coh(Z)
    Fz = (F5 - F5.mean(0)) / F5.std(0)
    with np.errstate(all='ignore'):
        rf = (C0.X0 @ Fz) / C0.M.sum(1, keepdims=True)          # approx corr of each site with each PC
    comm = (rf ** 2).sum(1)
    cov = x.notna().mean(1).values; sd = x.std(1).values
    res = SA.residue.reindex(Z.index).fillna('S').values; p1 = SA.pro1.reindex(Z.index).fillna(False).values
    q = lambda v, n: pd.qcut(pd.Series(v).rank(method='first'), n, labels=False).values
    binA = q(cov, 5) * 5 + q(sd, 5)
    binB = ((q(cov, 3) * 3 + q(sd, 3)) * 4 + q(comm, 4)) * 4 + (res == 'Y') * 2 + p1.astype(int)
    pools = {}
    for nm, b in [('A', binA), ('B', binB)]:
        pools[nm] = (b, {v: np.where(b == v)[0] for v in np.unique(b)})
    feats.append(pd.DataFrame({'site': Z.index, 'cohort': coh, 'coverage': cov, 'sd': sd, 'communality5': comm, 'residue': res, 'pro1': p1}))
    pidx = Pc.idx
    up = np.unique([s_.split('_')[0] for s_ in Z.index]); uix = {u: i for i, u in enumerate(up)}
    st = STR[STR.u1.isin(uix) & STR.u2.isin(uix)]
    ADJ = np.zeros((len(up), len(up)), bool); a_ = st.u1.map(uix).values; b2 = st.u2.map(uix).values
    ADJ[a_, b2] = True; ADJ[b2, a_] = True
    def extra(ii, kind, C):
        if kind is None: return None
        ps = C.prot[ii]; m = np.ones((len(ii), len(ii)), bool)
        if 'coab' in kind:
            pj = np.array([pidx.get(p, -1) for p in ps]); ok = pj >= 0
            rp = np.full((len(ii), len(ii)), np.nan); sub = np.where(ok)[0]
            rp[np.ix_(sub, sub)] = Pc.corr(pj[sub])
            m &= ~(np.abs(rp) > 0.3)
        if 'string' in kind:
            u = np.array([uix[p] for p in ps]); m &= ~ADJ[np.ix_(u, u)]
        return m
    conds = {'base': (C0, None), 'coab': (C0, 'coab'), 'string': (C0, 'string'), 'coab+string': (C0, 'coab+string'),
             'pc1': (Coh(Rk[1]), None), 'pc3': (Coh(Rk[3]), None), 'pc5': (Coh(Rk[5]), None), 'pc5+coab+string': (Coh(Rk[5]), 'coab+string')}
    for cname, (C, kind) in conds.items():
        for nm, sg in SETS.items():
            sg = sg[sg.index.isin(C.idx)]
            if len(sg) < 10: continue
            ii = np.array([C.idx[s] for s in sg.index]); sv = sg.values
            obs, npair = C.coherence(ii, sv, extra(ii, kind, C))
            row = dict(set=nm, cohort=coh, cond=cname, n_sites=len(sg), n_pairs=npair, coherence=obs,
                       cov_mean=cov[ii].mean(), sd_mean=sd[ii].mean(), comm_mean=comm[ii].mean(), frac_Y=(res[ii] == 'Y').mean())
            if nm in NULL_SETS and (cname in ('base', 'coab+string', 'pc5', 'pc5+coab+string') or nm.startswith('FP')):
                for pn, (b, pool) in pools.items():
                    nl = []
                    for _ in range(NN if cname in ('base', 'pc5+coab+string') else NN // 2):
                        jj = np.array([rng.choice(pool[v]) for v in b[ii]])
                        nl.append(C.coherence(jj, sv, extra(jj, kind, C))[0])
                    nl = np.array(nl)
                    row[f'null{pn}_mean'] = np.nanmean(nl); row[f'z{pn}'] = (obs - np.nanmean(nl)) / np.nanstd(nl)
                    row[f'p{pn}'] = (1 + (nl >= obs).sum()) / (1 + np.isfinite(nl).sum())
            rows.append(row)
    C5 = conds['pc5'][0]
    for e, sg in REXP.items():
        sg = sg[sg.index.isin(C0.idx)]
        if len(sg) < 30: continue
        ii = np.array([C0.idx[s] for s in sg.index])
        rrows.append(dict(exp=e, cohort=coh, n=len(sg), coh_base=C0.coherence(ii, sg.values)[0],
                          coh_pc5=C5.coherence(ii, sg.values)[0],
                          comm_mean=comm[ii].mean()))
    print(coh, round(time.time() - t0), 's', flush=True)
    pd.DataFrame(rows).to_csv(RA + 'A1_coherence.tsv', sep='\t', index=False)
pd.DataFrame(rrows).to_csv(RA + 'A1_random_experiments.tsv', sep='\t', index=False)
pd.concat(feats).to_parquet(CA + 'site_features.parquet')
