"""b_Q8 step 1: footprint coherence in NATs vs size-matched tumours, per cohort.
Groups: NAT = NATs of paired cases; T = tumours of the SAME cases (exact n and patient matching). UCEC (unpaired NX donors): T = random
tumour subsample of NAT n. Sensitivity: 5 random tumour subsamples (any patients) of the same n.
Conditions (all processing within group): raw | padj | padj_pc3 | padj_pc5 | padj_comp (stromal+immune+epithelial protein proxies +
ischaemia-footprint score regressed out) | padj_comp_pc3.
Statistics per set: signed mean cross-protein pairwise r (coherence); matched random-site null per group (bins: coverage x SD x
communality(top-5 PCs) x residue Y x proline+1; 100 draws, footprint signs kept) -> z_T, z_N;
differential Delta = C_T - C_N with (i) matched random-site Delta null (same random sites in both groups, bins from pooled features)
and (ii) label permutation (swap T/N within patient; 100 perms) -> p_perm, z_perm.
Outputs: tables/T1_coherence.tsv (+cache footprint scores for step 2)."""
import numpy as np, pandas as pd, warnings, time, sys
warnings.filterwarnings('ignore')
from common import *
rng = np.random.default_rng(8); NN = 100; NP = 100
SETS = {f'FP:{g}': fp_set(g) for g in GROUPS}
PS = psp_sets(15)
for k in ['CDK1', 'CDK2', 'MAPK1', 'MAPK3', 'CSNK2A1', 'PRKACA', 'AKT1', 'MTOR', 'ATM', 'PRKCA']:
    SETS[f'PSP:{k}'] = pd.Series(1.0, index=PS[k])
ISCH = ischaemia_set(); CS = pd.read_parquet(CQ + 'composition_scores.parquet')
SA = site_annot()
COND = {'raw': dict(padj=False, k=0), 'padj': dict(padj=True, k=0), 'padj_pc3': dict(padj=True, k=3), 'padj_pc5': dict(padj=True, k=5),
        'padj_comp': dict(padj=True, k=0, comp=True), 'padj_comp_pc3': dict(padj=True, k=3, comp=True)}
cohorts = sys.argv[1:] or NATC

def feats(Zraw_sd, Z, cov):
    F, _ = global_pcs(Z, 5); Fz = (F - F.mean(0)) / F.std(0)
    X0 = np.nan_to_num(Z.values); M = np.isfinite(Z.values)
    comm = (((X0 @ Fz) / np.maximum(M.sum(1, keepdims=True), 1)) ** 2).sum(1)
    res = SA.residue.reindex(Z.index).fillna('S').values; p1 = SA.pro1.reindex(Z.index).fillna(False).values.astype(int)
    q = lambda v, n: pd.qcut(pd.Series(v).rank(method='first'), n, labels=False).values
    b = ((q(cov, 3) * 3 + q(Zraw_sd, 3)) * 4 + q(comm, 4)) * 4 + (res == 'Y') * 2 + p1
    return b, {v: np.where(b == v)[0] for v in np.unique(b)}

def draw(b, pool, ii):
    return np.array([rng.choice(pool[v]) for v in b[ii]])

rows, scores = [], []
for coh in cohorts:
    t0 = time.time()
    ph, pr, T, N, paired = load_raw(coh)
    if len(paired) >= 15:
        Ns = [f'{c}_N' for c in paired if f'{c}_N' in N]; Ts = [f'{c}_T' for c in paired if f'{c}_T' in T]
        cs = sorted(set(c[:-2] for c in Ns) & set(c[:-2] for c in Ts)); Ns = [c + '_N' for c in cs]; Ts = [c + '_T' for c in cs]; design = 'paired'
    else:
        Ns = N; Ts = list(rng.choice(T, len(N), replace=False)); design = 'unpaired_subsample'
    n = len(Ns); minpair = max(10, int(0.5 * n))
    subs = [list(rng.choice(T, n, replace=False)) for _ in range(5)]
    print(coh, design, 'n', n, 'sites', len(ph), flush=True)
    for cname, cfg in COND.items():
        if cname == 'padj_pc5' and n < 40: continue
        def prep(samples):
            cov = None
            if cfg.get('comp'):
                Zr = standardise(ph[samples]); isc = score(Zr, ISCH)
                cov = CS.loc[samples, ['stromal', 'immune', 'epithelial']].assign(isch=isc).fillna(0)
                cov = cov.fillna(cov.mean())
            return process(ph, pr, samples, padj=cfg['padj'], k=cfg['k'], cov=cov, minobs=max(8, int(0.4 * n)))
        ZT, ZN = prep(Ts), prep(Ns)
        common_sites = ZT.index.intersection(ZN.index); ZT, ZN = ZT.loc[common_sites], ZN.loc[common_sites]
        CT, CN = Coh(ZT, minpair), Coh(ZN, minpair)
        cov_all = pd.concat([ph.loc[common_sites, Ts], ph.loc[common_sites, Ns]], axis=1)
        bT = feats(ph.loc[common_sites, Ts].std(1).values, ZT, ph.loc[common_sites, Ts].notna().mean(1).values)
        bN = feats(ph.loc[common_sites, Ns].std(1).values, ZN, ph.loc[common_sites, Ns].notna().mean(1).values)
        bP = feats(cov_all.std(1).values, pd.concat([ZT, ZN], axis=1), cov_all.notna().mean(1).values)
        # permutation matrices (swap within pairs)
        XT, XN = ZT.values, ZN.values
        perms = []
        if cname in ('padj_pc3', 'padj_comp_pc3', 'padj'):
            for _ in range(NP):
                if design == 'paired':
                    sw = rng.random(n) < 0.5
                    A = np.where(sw[None, :], XN, XT); B = np.where(sw[None, :], XT, XN)
                else:
                    P = np.concatenate([XT, XN], 1)[:, rng.permutation(2 * n)]; A, B = P[:, :n], P[:, n:]
                perms.append((Coh(pd.DataFrame(A, index=common_sites), minpair), Coh(pd.DataFrame(B, index=common_sites), minpair)))
        subC = []
        if cname == 'padj_pc3':
            for s_ in subs:
                Zs = prep(s_).reindex(common_sites); subC.append(Coh(Zs, minpair))
        for nm, sg in SETS.items():
            sg = sg[sg.index.isin(CT.idx)]
            if len(sg) < 15: continue
            ii = np.array([CT.idx[s] for s in sg.index]); sv = sg.values
            cT, npT = CT.coherence(ii, sv); cN, npN = CN.coherence(ii, sv)
            row = dict(set=nm, cohort=coh, cond=cname, design=design, n=n, n_sites=len(sg), n_pairs_T=npT, n_pairs_N=npN, coh_T=cT, coh_N=cN, delta=cT - cN)
            nT = np.array([CT.coherence(draw(*bT, ii), sv)[0] for _ in range(NN)])
            nN = np.array([CN.coherence(draw(*bN, ii), sv)[0] for _ in range(NN)])
            row.update(null_T=np.nanmean(nT), z_T=(cT - np.nanmean(nT)) / np.nanstd(nT), p_T=(1 + (nT >= cT).sum()) / (1 + NN),
                       null_N=np.nanmean(nN), z_N=(cN - np.nanmean(nN)) / np.nanstd(nN), p_N=(1 + (nN >= cN).sum()) / (1 + NN))
            dn = []
            for _ in range(NN):
                jj = draw(*bP, ii); dn.append(CT.coherence(jj, sv)[0] - CN.coherence(jj, sv)[0])
            dn = np.array(dn); row.update(null_delta=np.nanmean(dn), z_delta_null=(row['delta'] - np.nanmean(dn)) / np.nanstd(dn))
            if perms:
                pd_ = np.array([a.coherence(ii, sv)[0] - b.coherence(ii, sv)[0] for a, b in perms])
                row.update(z_delta_perm=(row['delta'] - np.nanmean(pd_)) / np.nanstd(pd_),
                           p_delta_perm=(1 + (np.abs(pd_ - np.nanmean(pd_)) >= abs(row['delta'] - np.nanmean(pd_))).sum()) / (1 + NP))
            if subC:
                row['coh_Tsub_mean'] = np.nanmean([c.coherence(np.array([c.idx[s] for s in sg.index]), sv)[0] for c in subC])
            rows.append(row)
        if cname in ('padj', 'padj_pc3', 'padj_comp_pc3'):
            for grp, Z in [('T', ZT), ('N', ZN)]:
                sc = pd.DataFrame({nm: score(Z, sg) for nm, sg in SETS.items()})
                sc['cohort'] = coh; sc['grp'] = grp; sc['cond'] = cname; scores.append(sc)
        print(' ', cname, round(time.time() - t0), 's', flush=True)
    pd.DataFrame(rows).to_csv(TQ + f'T1_coherence_{"_".join(cohorts) if len(cohorts) < 8 else "all"}.tsv', sep='\t', index=False)
    pd.concat(scores).to_parquet(CQ + f'fp_scores_{"_".join(cohorts) if len(cohorts) < 8 else "all"}.parquet')
