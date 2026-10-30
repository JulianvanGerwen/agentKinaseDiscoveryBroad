"""I10 step 2: consensus site-level perturbation footprints + held-out-study validation.
Unit = study (groups with >=2 studies) else cell line. Unit vector = mean of its experiments' robust z.
Site kept if measured in >=2 units, sign-consistent in >=75% of units (100% if 2 units), |mean z|>=1;
weight = Stouffer S = sum z_u / sqrt(n_u) (clipped at +-8); footprint = top K by |S| (K=200), plus down-only / up-only variants.
Validation: leave-one-study-out (LOSO): footprints rebuilt without study s, every labelled experiment of s scored
(cosine on overlapping sites); AUROC same-group vs other-group experiments (and vs other-axis only)."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
K = 200
M = pd.read_csv(CA + 'fp_meta.tsv', sep='\t', dtype={'study': str}).set_index('exp')
W = pd.read_parquet(CA + 'fp_sites.parquet')[M.index]
G = pd.read_csv(RT + 'T0_groups.tsv', sep='\t').set_index('group')
elig = G.index[G.eligible].tolist()
# composite groups
comp = {'MAPKaxis': ['KI_MEK', 'KI_ERK', 'KI_RAF'], 'PI3Kaxis': ['KI_PI3K', 'KI_PI3K_mTOR', 'KI_mTORcat', 'KI_AKT', 'KI_PDPK1', 'KI_rapalog', 'KI_S6K'],
        'GF_all': ['EGF', 'FGF', 'insulin', 'serum', 'starvation_rev']}
M['sgn'] = 1.0
st = M[M.group == 'starvation'].copy(); st['group'] = 'starvation_rev'; st['sgn'] = -1.0; st.index = st.index + '#rev'
W = pd.concat([W, -W[[i[:-4] for i in st.index]].set_axis(st.index, axis=1)], axis=1)
M = pd.concat([M, st])
members = {g: M.index[M.group == g].tolist() for g in elig}
for c, gs in comp.items(): members[c] = M.index[M.group.isin(gs)].tolist()
Wv = W.values.astype(np.float32); Wv = np.clip(Wv, -10, 10)
col = {e: i for i, e in enumerate(W.columns)}
sites = W.index.values

def unit_matrix(exps):
    m = M.loc[exps]; byst = m.study.nunique() >= 2
    key = m.study if byst else m.cell.astype(str)
    U, ids = [], []
    for u, d in m.groupby(key):
        U.append(np.nanmean(Wv[:, [col[e] for e in d.index]], axis=1)); ids.append(u)
    return np.column_stack(U), ids, byst

def build(U):
    n = np.isfinite(U).sum(1)
    with np.errstate(all='ignore'):
        mean = np.nanmean(U, 1); S = np.nansum(U, 1) / np.sqrt(n)
        cons = (np.sign(U) == np.sign(mean)[:, None]).sum(1) / n
    nu = U.shape[1]
    ok = (n >= 2) & (np.abs(mean) >= 1) & (cons >= (1.0 if nu == 2 else 0.75))
    if nu == 1: ok = np.isfinite(U[:, 0]) & (np.abs(U[:, 0]) >= 2); S = np.nan_to_num(U[:, 0]); mean = S; cons = np.ones_like(S); n = np.isfinite(U[:, 0]).astype(int)
    idx = np.where(ok)[0]
    S = np.clip(S, -8, 8)
    top = idx[np.argsort(-np.abs(S[idx]))][:K]
    return pd.DataFrame({'site': sites[top], 'w': S[top], 'mean_z': mean[top], 'n_units': n[top], 'consistency': cons[top]})

def score(fp, exps, which='all'):
    f = fp if which == 'all' else fp[fp.w < 0] if which == 'down' else fp[fp.w > 0]
    if len(f) < 5: return pd.Series(np.nan, index=exps)
    r = W.index.get_indexer(f.site); w = f.w.values
    out = {}
    for e in exps:
        z = Wv[r, col[e]]; m = np.isfinite(z)
        if m.sum() < 10: out[e] = np.nan; continue
        z = np.clip(z[m], -5, 5); ww = w[m]
        if which == 'all': out[e] = (ww * z).sum() / np.sqrt((ww ** 2).sum() * (z ** 2).sum() + 1e-9)
        else: out[e] = (np.abs(ww) * z * np.sign(ww)).sum() / np.abs(ww).sum()   # signed weighted mean (direction-specific)
    return pd.Series(out)

if __name__ == '__main__':
    # final footprints
    FP = []
    for g, ex in members.items():
        U, ids, byst = unit_matrix(ex); f = build(U); f['group'] = g; f['unit'] = 'study' if byst else 'cell'; f['n_units_total'] = len(ids); FP.append(f)
    FP = pd.concat(FP); FP.to_csv(RT + 'T1_footprints.tsv', sep='\t', index=False)
    s = FP.groupby('group').agg(n_sites=('site', 'size'), n_down=('w', lambda x: (x < 0).sum()), unit=('unit', 'first'), n_units=('n_units_total', 'first'))
    print(s.to_string())
    # LOSO
    multi = [g for g in members if M.loc[members[g]].study.nunique() >= 2]
    lab = M[M.group.isin(elig) & ~M.index.str.endswith('#rev')]
    rows = []
    for s_ in sorted(lab.study.unique()):
        test = lab.index[lab.study == s_].tolist()
        for g in multi:
            tr = [e for e in members[g] if M.loc[e, 'study'] != s_]
            if not tr: continue
            U, ids, byst = unit_matrix(tr); f = build(U)
            for which in ['all', 'down', 'up']:
                sc = score(f, test, which)
                for e, v in sc.items(): rows.append((s_, e, M.loc[e, 'group'], g, which, v, len(ids)))
        print('LOSO', s_, len(test), flush=True)
    R = pd.DataFrame(rows, columns=['heldout_study', 'exp', 'exp_group', 'footprint', 'variant', 'score', 'n_train_units'])
    R.to_parquet(CA + 'loso_scores.parquet')
    # AUROC per footprint
    ax = lambda g: next((c for c, gs in comp.items() if g in gs or g == c), g)
    out = []
    for (fpg, v), d in R.dropna(subset=['score']).groupby(['footprint', 'variant']):
        pos_groups = comp.get(fpg, [fpg])
        y = d.exp_group.isin(pos_groups)
        if y.sum() < 2 or (~y).sum() < 5: continue
        d2 = d[y | (d.exp_group.map(ax) != ax(fpg))]; y2 = d2.exp_group.isin(pos_groups)
        out.append((fpg, v, int(y.sum()), int((~y).sum()), roc_auc_score(y, d.score), roc_auc_score(y2, d2.score), d[y].heldout_study.nunique()))
    A = pd.DataFrame(out, columns=['footprint', 'variant', 'n_pos', 'n_neg', 'auroc', 'auroc_vs_other_axis', 'n_heldout_studies'])
    A = A.pivot(index='footprint', columns='variant', values=['auroc', 'auroc_vs_other_axis']).join(A[A.variant == 'all'].set_index('footprint')[['n_pos', 'n_neg', 'n_heldout_studies']])
    A.columns = ['_'.join(c) if isinstance(c, tuple) else c for c in A.columns]
    A.to_csv(RT + 'T2_loso_auroc.tsv', sep='\t'); print(A.round(3).to_string())
