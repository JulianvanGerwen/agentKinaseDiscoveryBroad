"""b_Q8 step 4: preservation of compendium-derived relationships in NATs vs tumours (same paired design as step 1).
(a) programme->programme pairs (T3: footprints + PSP kinase sets; exclusive site sets): r_T, r_N per cohort.
(b) programme->site edges: each footprint member site vs the footprint score computed without any site on that protein (sign = footprint w).
Differential: Fisher z_diff = (atanh r_T - atanh r_N)/sqrt(1/(nT-3)+1/(nN-3)); permutation (swap T/N within patient, 200 perms; scores are
per-sample linear functionals of within-group processed data, so permuting score vectors == permuting samples) -> p_perm, z_perm.
Conditions: padj_pc3 (primary), padj, padj_comp_pc3. Output: cache edge_stats_<cond>.parquet"""
import numpy as np, pandas as pd, warnings, time, sys
warnings.filterwarnings('ignore')
from common import *
from pairs import prepare_sets, pair_scores
rng = np.random.default_rng(5); NP = 200
E = pd.read_csv(TQ + 'T3_compendium_edges.tsv', sep='\t')
sets = {f'FP:{g}': fp_set(g) for g in GROUPS}
for k, v in psp_sets(15).items(): sets[f'K:{k}'] = pd.Series(1.0, index=v)
ISCH = ischaemia_set(); CS = pd.read_parquet(CQ + 'composition_scores.parquet')
COND = {'padj_pc3': dict(padj=True, k=3), 'padj': dict(padj=True, k=0), 'padj_comp_pc3': dict(padj=True, k=3, comp=True)}
conds = sys.argv[1:] or list(COND)

def corr_cols(A, B):
    """A, B: (n_samples x m) arrays with NaN; column-wise pairwise-complete Pearson."""
    M = np.isfinite(A) & np.isfinite(B); n = M.sum(0)
    A = np.where(M, A, np.nan); B = np.where(M, B, np.nan)
    with np.errstate(all='ignore'):
        A = A - np.nanmean(A, 0); B = B - np.nanmean(B, 0)
        r = np.nansum(A * B, 0) / np.sqrt(np.nansum(A ** 2, 0) * np.nansum(B ** 2, 0))
    r[n < 12] = np.nan; return r, n

for cname in conds:
    cfg = COND[cname]; out = []; gl = []
    for coh in NATC:
        t0 = time.time()
        ph, pr, T, N, paired = load_raw(coh)
        if len(paired) >= 15:
            cs = sorted(set(c for c in paired if f'{c}_N' in N and f'{c}_T' in T)); Ns = [c + '_N' for c in cs]; Ts = [c + '_T' for c in cs]; design = 'paired'
        else:
            Ns = N; Ts = list(np.random.default_rng(8).choice(T, len(N), replace=False)); design = 'unpaired'
        n = len(Ns)
        def prep(samples):
            cov = None
            if cfg.get('comp'):
                isc = score(standardise(ph[samples]), ISCH)
                cov = CS.loc[samples, ['stromal', 'immune', 'epithelial']].assign(isch=isc); cov = cov.fillna(cov.mean())
            return process(ph, pr, samples, padj=cfg['padj'], k=cfg['k'], cov=cov, minobs=max(8, int(0.4 * n)))
        ZT, ZN = prep(Ts), prep(Ns)
        cs_ = ZT.index.intersection(ZN.index); Z = pd.concat([ZT.loc[cs_], ZN.loc[cs_]], axis=1)   # T cols then N cols
        S = prepare_sets(sets, Z.index)
        # (a) programme pairs
        pr_ = [(a, b) for a, b in zip(E.a, E.b) if a in S and b in S and len(S[a]) >= 10 and len(S[b]) >= 10]
        PS = pair_scores(Z, S, pr_)
        keys = [k for k, v in PS.items() if v[0] is not None and v[1] is not None]
        A = np.column_stack([PS[k][0] for k in keys]); B = np.column_stack([PS[k][1] for k in keys])
        lab = [f'{a}|{b}' for a, b in keys]; kind = ['pair'] * len(keys)
        # (b) footprint -> member site
        X0 = np.nan_to_num(Z.values); Mk = np.isfinite(Z.values).astype(float); idx = {s: i for i, s in enumerate(Z.index)}
        As, Bs, ls = [], [], []
        for g in GROUPS:
            s = S[f'FP:{g}']
            if len(s) < 20: continue
            ii = np.array([idx[x] for x in s.index]); w = s.values; prot = np.array([x.split('_')[0] for x in s.index])
            num = (X0[ii] * w[:, None]); cnt = Mk[ii]
            tot, totn = num.sum(0), cnt.sum(0)
            for p in np.unique(prot):
                jj = prot == p; sc = (tot - num[jj].sum(0)) / np.maximum(totn - cnt[jj].sum(0), 1)
                sc[(totn - cnt[jj].sum(0)) < 10] = np.nan
                for k in np.where(jj)[0]:
                    As.append(Z.values[ii[k]] * w[k]); Bs.append(sc); ls.append(f'FP:{g}>{s.index[k]}')
        A = np.column_stack([A] + [np.array(As).T]); B = np.column_stack([B] + [np.array(Bs).T])
        lab += ls; kind += ['site'] * len(ls)
        At, Bt, An, Bn = A[:n], B[:n], A[n:], B[n:]
        rT, nT = corr_cols(At, Bt); rN, nN = corr_cols(An, Bn)
        zd = (np.arctanh(np.clip(rT, -.999, .999)) - np.arctanh(np.clip(rN, -.999, .999))) / np.sqrt(1 / (nT - 3) + 1 / (nN - 3))
        emap = {f'{a}|{b}': np.sign(r) for a, b, r, l in zip(E.a, E.b, E.r_meta, E.label) if l == 'edge'}
        sgn = np.array([emap.get(l, 0.0) for l in lab]); isk = np.array(kind) == 'site'
        def glob(rx, ry):
            dd = rx - ry
            return [np.nanmean((sgn * dd)[sgn != 0]), np.nanmean(dd[isk]), np.nanmean((sgn * rx)[sgn != 0]), np.nanmean((sgn * ry)[sgn != 0]), np.nanmean(rx[isk]), np.nanmean(ry[isk])]
        G = [['obs'] + glob(rT, rN)]
        dP = []
        for _ in range(NP):
            if design == 'paired':
                sw = rng.random(n) < 0.5; sel1 = np.where(sw, np.arange(n) + n, np.arange(n)); sel2 = np.where(sw, np.arange(n), np.arange(n) + n)
            else:
                p_ = rng.permutation(2 * n); sel1, sel2 = p_[:n], p_[n:]
            r1, _ = corr_cols(A[sel1], B[sel1]); r2, _ = corr_cols(A[sel2], B[sel2]); dP.append(r1 - r2); G.append(['perm'] + glob(r1, r2))
        dP = np.array(dP); d = rT - rN
        mu, sd = np.nanmean(dP, 0), np.nanstd(dP, 0)
        pperm = (1 + (np.abs(dP - mu) >= np.abs(d - mu)).sum(0)) / (1 + NP)
        out.append(pd.DataFrame(dict(rel=lab, kind=kind, cohort=coh, n=n, r_T=rT, r_N=rN, nT=nT, nN=nN, z_fisher=zd, z_perm=(d - mu) / sd, p_perm=pperm)))
        gl.append(pd.DataFrame(G, columns=['what', 'dpres_pairs', 'dpres_sites', 'pres_pairs_T', 'pres_pairs_N', 'pres_sites_T', 'pres_sites_N']).assign(cohort=coh))
        print(cname, coh, len(lab), round(time.time() - t0), 's', flush=True)
    pd.concat(out).to_parquet(CQ + f'edge_stats_{cname}.parquet'); pd.concat(gl).to_parquet(CQ + f'edge_global_{cname}.parquet')
