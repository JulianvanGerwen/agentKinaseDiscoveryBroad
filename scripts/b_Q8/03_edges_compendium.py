"""b_Q8 step 3: compendium-derived signed programme->programme relationships (Invergo-style co-regulation across perturbations).
Nodes: 49 consensus perturbation footprints (signed site sets, b_I10) + PSP kinase substrate sets (>=15 sites measured in CPTAC union).
Relationship: for each pair, exclusive site sets (no shared sites, no sites on proteins carrying the other set's sites); per-experiment
scores = signed mean z (688 b_I10 consensus experiments, >=5 sites). Per study (>=6 experiments with both scored) Pearson r,
Fisher-z meta (weights n-3, capped at 40 per study to stop decryptM dominating), sign agreement across studies.
Edge = |meta r| >= 0.3, >= 3 studies, >= 70% studies sign-concordant, p_meta < 1e-4. Non-edge baseline: |meta r| < 0.05, >= 3 studies.
Output: tables/T3_compendium_edges.tsv (all tested pairs with label)."""
import numpy as np, pandas as pd, itertools, warnings
warnings.filterwarnings('ignore')
from common import *
from pairs import prepare_sets, pair_scores, rcorr
from scipy.stats import norm
W = pd.read_parquet('/workspace/cache/b_I10/fp_sites.parquet')
# universe = sites measured in >=1 CPTAC NAT cohort site universe (keeps nodes comparable to tumour tests)
univ = set()
for coh in MAINC:
    ph = cl.matrix(coh, 'phospho'); T, N, _ = groups(coh)
    ok = lambda X: (X.notna().mean(1) >= 0.25) & (X.notna().sum(1) >= 10)
    univ |= set(ph.index[ok(ph[[c for c in T if c in ph.columns]]) & ok(ph[[c for c in N if c in ph.columns]])])
W = W[W.index.isin(univ)]
print('compendium x CPTAC sites', W.shape)
sets = {f'FP:{g}': fp_set(g) for g in GROUPS}
for k, v in psp_sets(15).items(): sets[f'K:{k}'] = pd.Series(1.0, index=v)
sets = prepare_sets(sets, W.index)
sets = {k: v for k, v in sets.items() if len(v) >= 15}
print('nodes', len(sets), sum(k.startswith('K:') for k in sets))
# robust-standardise each experiment
Z = (W - W.median()) / (1.4826 * (W - W.median()).abs().median())
study = np.array([c.split('|')[0] for c in Z.columns])
pairs = list(itertools.combinations(sorted(sets), 2))
PS = pair_scores(Z, sets, pairs)
rows = []
for (a, b), (sa, sb) in PS.items():
    if sa is None or sb is None: continue
    zs, ws, sgn, nst = [], [], [], 0; ntot = 0
    for st in np.unique(study):
        m = study == st
        r, nn = rcorr(sa[m], sb[m], 6)
        if not np.isfinite(r) or nn < 6: continue
        zs.append(np.arctanh(np.clip(r, -0.999, 0.999))); ws.append(min(nn, 40) - 3); sgn.append(np.sign(r)); ntot += nn
    if len(zs) < 3: continue
    zs, ws = np.array(zs), np.array(ws)
    zm = (zs * ws).sum() / ws.sum(); se = 1 / np.sqrt(ws.sum())
    rp, _ = rcorr(sa, sb)
    rows.append(dict(a=a, b=b, n_studies=len(zs), n_exp=ntot, r_meta=np.tanh(zm), z_meta=zm / se, p_meta=2 * norm.sf(abs(zm / se)),
                     frac_sign=np.mean(np.array(sgn) == np.sign(zm)), r_pooled=rp))
E = pd.DataFrame(rows)
E['label'] = 'other'
E.loc[(E.r_meta.abs() >= 0.3) & (E.n_studies >= 3) & (E.frac_sign >= 0.7) & (E.p_meta < 1e-4), 'label'] = 'edge'
E.loc[(E.r_meta.abs() < 0.05) & (E.n_studies >= 3), 'label'] = 'nonedge'
E['type'] = [('FP' if a.startswith('FP') else 'K') + '-' + ('FP' if b.startswith('FP') else 'K') for a, b in zip(E.a, E.b)]
E.to_csv(TQ + 'T3_compendium_edges.tsv', sep='\t', index=False)
print(E.groupby(['type', 'label']).size().unstack())
print(E[E.label == 'edge'].sort_values('z_meta').head(10)[['a', 'b', 'r_meta', 'n_studies', 'frac_sign']])
print(E[E.label == 'edge'].sort_values('z_meta').tail(10)[['a', 'b', 'r_meta', 'n_studies', 'frac_sign']])
