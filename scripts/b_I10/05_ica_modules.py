"""I10 step 5: ICA phosphomodules on the compendium (iModulon-style).
Matrix: module-experiments (robust z, clipped +-6) x sites measured in >=30% of experiments; experiments with >=30% of these sites;
missing -> 0 (no change vs control). FastICA (sites = samples, so sources are site-weight vectors), k=40,
20 random restarts on the full matrix + 20 bootstrap resamples of experiments. Components clustered across runs
(|Pearson| of site weights >= 0.7, sign-aligned); robust = cluster present in >=50% of runs. Module members: |robust z of weight| >= 4 (max 100).
Outputs: cache/b_I10/ica_runs.npz, results tables T6_modules.tsv, T6_module_sites.tsv, cache/b_I10/module_activity.parquet"""
import numpy as np, pandas as pd, warnings
from sklearn.decomposition import FastICA
warnings.filterwarnings('ignore')
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
import os
K = 40; NRUN = 20
CENTER = os.environ.get('ICA_CENTER', '0') == '1'; SUF = '_sc' if CENTER else ''
W = pd.read_parquet(CA + 'mod_sites.parquet')
M = pd.read_csv(CA + 'mod_meta.tsv', sep='\t', dtype={'study': str}).set_index('exp')
W = W[W.notna().mean(1) >= 0.30]
W = W.loc[:, W.notna().mean(0) >= 0.30]
if CENTER:  # within-study centring (studies with >=5 experiments): removes study/plex-shared shifts
    st = M.loc[W.columns, 'study']
    for s_, cols in st.groupby(st).groups.items():
        if len(cols) >= 5: W[cols] = W[cols].sub(W[cols].median(1), axis=0)
print('matrix', W.shape, 'studies', M.loc[W.columns].study.nunique())
X = np.clip(W.fillna(0).values, -6, 6).astype(np.float64)
X = X - X.mean(0)  # centre each experiment
comps = []; runlab = []
for r in range(2 * NRUN):
    rng = np.random.default_rng(r)
    cols = np.arange(X.shape[1]) if r < NRUN else rng.choice(X.shape[1], X.shape[1], replace=True)
    ica = FastICA(n_components=K, random_state=r, max_iter=1000, tol=1e-4, whiten='unit-variance')
    S = ica.fit_transform(X[:, cols])            # sites x K
    S = S * np.sign(pd.DataFrame(S).skew().values)  # orient heavy tail positive
    comps.append(S); runlab += [r] * K
    print('run', r, flush=True)
S_all = np.hstack(comps); runlab = np.array(runlab)
np.savez_compressed(CA + 'ica_runs' + SUF + '.npz', S=S_all.astype(np.float32), run=runlab, sites=W.index.values)
C = np.corrcoef(S_all.T); A = np.abs(C)
# greedy clustering: seeds = components with most partners
n = A.shape[0]; assigned = -np.ones(n, int); clusters = []
deg = (A >= 0.7).sum(1)
for i in np.argsort(-deg):
    if assigned[i] >= 0: continue
    mem = np.where((A[i] >= 0.7) & (assigned < 0))[0]
    # one component per run
    keep = []; seen = set()
    for j in mem[np.argsort(-A[i, mem])]:
        if runlab[j] not in seen: keep.append(j); seen.add(runlab[j])
    assigned[keep] = len(clusters); clusters.append(keep)
rob = [c for c in clusters if len(c) >= NRUN]   # >=50% of 40 runs
print('clusters', len(clusters), 'robust', len(rob))
mods, msites = [], []
for m, c in enumerate(rob):
    ref = S_all[:, c[0]]
    V = np.column_stack([S_all[:, j] * np.sign(np.corrcoef(ref, S_all[:, j])[0, 1]) for j in c])
    w = V.mean(1)
    z = (w - np.median(w)) / (np.median(np.abs(w - np.median(w))) * 1.4826)
    if (z > 0).sum() and abs(z.min()) > z.max(): z = -z      # dominant tail positive
    order = np.argsort(-np.abs(z)); mem = [i for i in order[:100] if abs(z[i]) >= 4]
    stab = np.mean(A[np.ix_(c, c)][np.triu_indices(len(c), 1)])
    mods.append(dict(module=f'M{m+1:02d}', n_runs=len(c), stability=stab, n_sites=len(mem), n_pos=int((z[mem] > 0).sum())))
    for i in mem: msites.append((f'M{m+1:02d}', W.index[i], z[i]))
    if m == 0: Z = np.zeros((len(W), len(rob)))
    Z[:, m] = z
mods = pd.DataFrame(mods); MS = pd.DataFrame(msites, columns=['module', 'site', 'weight_z'])
# module activity per experiment: least squares of experiment vectors on module weight vectors (members only)
Zm = pd.DataFrame(Z, index=W.index, columns=mods.module)
Zm.astype('float32').to_parquet(CA + 'module_weights' + SUF + '.parquet')
act = {}
for mname, d in MS.groupby('module'):
    w = d.set_index('site').weight_z; sub = W.loc[w.index]
    act[mname] = (sub.mul(np.sign(w), axis=0)).mean(0)  # signed mean z of members per experiment (measured only)
act = pd.DataFrame(act); act.astype('float32').to_parquet(CA + 'module_activity' + SUF + '.parquet')
# which perturbations drive each module + study concentration
desc = []
for mname in mods.module:
    a = act[mname].dropna(); top = a.abs().sort_values(ascending=False).index[:20]
    mm = M.loc[top]
    st = mm.study.value_counts(); gr = mm.group.fillna(mm.subclass).value_counts()
    desc.append(dict(module=mname, top_study=st.index[0], top_study_frac=st.iloc[0] / 20, n_studies_top20=len(st),
                     top_groups=';'.join(f'{k}:{v}' for k, v in gr.head(4).items()),
                     top_treats=';'.join(mm.treat.astype(str).head(6))))
mods = mods.merge(pd.DataFrame(desc), on='module')
mods.to_csv(RT + 'T6_modules' + SUF + '.tsv', sep='\t', index=False); MS.to_csv(RT + 'T6_module_sites' + SUF + '.tsv', sep='\t', index=False)
print(mods.to_string())
