"""Task 2 + 4: kinase-level compensation atlas.
For each inhibitor class K and each kinase X (n_substrates>=5), estimate the change in X's inferred activity under K,
excluding experiments where X is a primary/Klaeger(<=1uM)/Hijazi target of the drug.
Values: within study x cell-line centred z (groups with >=5 exps and >=3 classes), otherwise vs DMSO.
Per-study mean -> DerSimonian-Laird random effects across studies (study = random effect).
Nulls: (a) matched random-condition null (same study counts, other-class experiments); (b) drug-identity shuffle within study x cell line.
LOO: leave-one-study-out. Replication: corev3 z_core vs psp_z (independent library)."""
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import *
OUT = "/workspace/results/b_I19/"
rng = np.random.default_rng(1)
E = experiments(); exps = E.exp.values
MIN_ST = 2
def load_act(name):
    M = pd.read_parquet(CA + f'act_{name}_exp.parquet').reindex(columns=exps)
    M = M[M.notna().sum(1) >= 20]
    Mc = center(M, E)
    B = np.array([[k in b for b in E.bset] for k in Mc.index])  # binder mask
    A = Mc.values.astype(float); A[B] = np.nan
    s = np.nanmedian(np.abs(Mc.values - np.nanmedian(Mc.values, 1, keepdims=True)), 1) * 1.4826
    return Mc.index.values, A, np.maximum(s, 0.3)
studies = E.study.values; ust = np.unique(studies)
def meta_vec(A, s, idx):
    """A: kin x exp; idx: boolean exp mask for class. returns DataFrame-like arrays per kinase"""
    sub = A[:, idx]; st = studies[idx]
    us = np.unique(st)
    Y = np.full((A.shape[0], len(us)), np.nan); N = np.zeros_like(Y)
    for j, u in enumerate(us):
        m = st == u; v = sub[:, m]; n = np.sum(~np.isnan(v), 1)
        with np.errstate(all='ignore'): Y[:, j] = np.where(n > 0, np.nanmean(v, 1), np.nan)
        N[:, j] = n
    V = (s[:, None] ** 2) / np.where(N > 0, N, np.nan)
    W = 1 / V; k = np.sum(~np.isnan(Y), 1)
    with np.errstate(all='ignore'):
        fe = np.nansum(W * Y, 1) / np.nansum(W, 1)
        Q = np.nansum(W * (Y - fe[:, None]) ** 2, 1); c = np.nansum(W, 1) - np.nansum(W ** 2, 1) / np.nansum(W, 1)
        tau2 = np.where(k > 1, np.maximum(0, (Q - (k - 1)) / c), 0)
        Ws = 1 / (V + tau2[:, None]); est = np.nansum(Ws * Y, 1) / np.nansum(Ws, 1); se = np.sqrt(1 / np.nansum(Ws, 1))
        I2 = np.where(Q > 0, np.maximum(0, (Q - (k - 1)) / Q), 0)
    fu = np.nanmean(np.where(np.isnan(sub), np.nan, sub > 0), 1)
    return dict(est=est, se=se, z=est / se, k=k, n=np.sum(~np.isnan(sub), 1), tau2=tau2, I2=I2, frac_up=fu, Y=Y, us=us)
classes = [c for c, d in E.groupby('cls') if d.study.nunique() >= MIN_ST]
res = {}
for meth in ['corev3', 'psp_z']:
    kin, A, s = load_act(meth)
    rows = []
    NP = 300
    for K in classes:
        idx = (E.cls == K).values
        r = meta_vec(A, s, idx)
        # null (a): matched random-condition draws from other classes in same study (fallback: global pool)
        nullA = np.zeros((NP, len(kin)))
        other = ~idx
        for p in range(NP):
            pick = np.zeros(len(E), bool)
            for u in np.unique(studies[idx]):
                nK = np.sum(idx & (studies == u)); pool = np.where(other & (studies == u))[0]
                if len(pool) < nK: pool = np.where(other)[0]
                pick[rng.choice(pool, nK, replace=False)] = True
            # keep study labels of the sampled experiments' *target* study for structure: use their own studies
            nullA[p] = meta_vec(A, s, pick)['est']
        # null (b): shuffle class labels within study x cell line groups
        nullB = np.zeros((NP, len(kin)))
        for p in range(NP):
            lab = E.cls.values.copy()
            for g, ii in E.groupby('grp').indices.items(): lab[ii] = rng.permutation(lab[ii])
            nullB[p] = meta_vec(A, s, lab == K)['est']
        with np.errstate(all='ignore'):
            pA = (1 + np.sum(nullA >= r['est'][None], 0)) / (NP + 1); pA_dn = (1 + np.sum(nullA <= r['est'][None], 0)) / (NP + 1)
            pB = (1 + np.sum(nullB >= r['est'][None], 0)) / (NP + 1); pB_dn = (1 + np.sum(nullB <= r['est'][None], 0)) / (NP + 1)
        # LOO
        loo_min = np.full(len(kin), np.inf); loo_max = np.full(len(kin), -np.inf)
        for u in np.unique(studies[idx]):
            rr = meta_vec(A, s, idx & (studies != u))
            e = np.where(np.isnan(rr['est']), r['est'], rr['est'])
            loo_min = np.minimum(loo_min, e); loo_max = np.maximum(loo_max, e)
        for i, X in enumerate(kin):
            if r['k'][i] == 0: continue
            rows.append(dict(cls=K, kinase=X, method=meth, est=r['est'][i], se=r['se'][i], z=r['z'][i], k_studies=int(r['k'][i]), n_exp=int(r['n'][i]),
                             frac_up=r['frac_up'][i], tau2=r['tau2'][i], I2=r['I2'][i], p_null_match_up=pA[i], p_null_match_dn=pA_dn[i],
                             p_null_shuffle_up=pB[i], p_null_shuffle_dn=pB_dn[i], loo_min=loo_min[i], loo_max=loo_max[i],
                             per_study=';'.join(f'{u}:{y:.2f}' for u, y in zip(r['us'], r['Y'][i]) if not np.isnan(y))))
        print(meth, K, 'done', flush=True)
    R = pd.DataFrame(rows)
    R['p'] = 2 * stats.norm.sf(R.z.abs())
    R['q'] = stats.false_discovery_control(R.p.fillna(1))
    res[meth] = R
R = pd.concat(res.values())
R.to_csv(OUT + 'tables/T2_compensation_atlas_all.tsv', sep='\t', index=False)
