"""I10 step 7: do compendium-derived site sets (ICA modules, footprints) form coherent co-varying units in CPTAC tumours?
Per cohort, protein-adjusted tumour matrix (tc_padj rank-normal, cached by step 3); sites centred/scaled, NaN->0 (mean imputation).
Signed coherence = mean over cross-protein pairs of s_i s_j r_ij (s = compendium sign). Sets need >=8 sites in the cohort.
Nulls: (a) 300 random site sets matched for (coverage x SD) bins with the same sign pattern; (b) the same sites with shuffled signs (300x).
Positive control: PSP kinase-substrate sets (>=10 sites in cohort), signs all +.
Output: T8_coherence.tsv (set x cohort: coh, z vs matched null, p_matched, p_signshuffle)."""
import sys, os, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
NN = 300; rng = np.random.default_rng(7)
SETS = {}
for suf, tag in [('', 'ICA'), ('_sc', 'ICAsc')]:
    if not os.path.exists(RT + f'T6_module_sites{suf}.tsv'): continue
    MS = pd.read_csv(RT + f'T6_module_sites{suf}.tsv', sep='\t')
    for m, d in MS.groupby('module'): SETS[f'{tag}:{m}'] = pd.Series(np.sign(d.weight_z.values), index=d.site.values)
FP = pd.read_csv(RT + 'T1_footprints.tsv', sep='\t')
for g, d in FP.groupby('group'):
    SETS[f'FP:{g}'] = pd.Series(np.sign(d.w.values), index=d.site.values)
S = load.sites(); ks = {}
for s, v in S.set_index('site').psp_kinases_human.dropna().items():
    for k in v.split(';'): ks.setdefault(k, []).append(s)
for k, v in ks.items():
    if len(v) >= 15: SETS[f'PSP:{k}'] = pd.Series(1.0, index=v)

def coh_score(Zc, prot, ii, sg):
    C = Zc[ii] @ Zc[ii].T
    P = np.outer(sg, sg) * C
    mask = np.triu(np.ones_like(C, bool), 1) & (prot[ii][:, None] != prot[ii][None, :])
    return P[mask].mean() if mask.sum() >= 10 else np.nan

rows = []
for coh in cl.COHORTS:
    R = pd.read_parquet(CA + f'cptac_rank_tc_padj_{coh}.parquet')
    X = R.values.astype(np.float64)
    mu = np.nanmean(X, 1, keepdims=True); sd = np.nanstd(X, 1, keepdims=True)
    Z = np.nan_to_num((X - mu) / sd); n = X.shape[1]
    Zc = Z / np.sqrt(n)
    cov = np.isfinite(X).mean(1); sdr = pd.Series(np.nanstd(X, 1))
    b = (pd.qcut(pd.Series(cov).rank(method='first'), 5, labels=False) * 5 + pd.qcut(sdr.rank(method='first'), 5, labels=False)).values
    bins = {k: np.where(b == k)[0] for k in range(25)}
    prot = np.array([s.split('_')[0] for s in R.index]); idx = {s: i for i, s in enumerate(R.index)}
    for nm, sg in SETS.items():
        sg = sg[sg.index.isin(idx)]
        if nm.startswith('PSP:') and len(sg) < 10: continue
        if len(sg) < 8: continue
        ii = np.array([idx[s] for s in sg.index]); sv = sg.values
        obs = coh_score(Zc, prot, ii, sv)
        nm_ = np.array([coh_score(Zc, prot, np.array([rng.choice(bins[q]) for q in b[ii]]), sv) for _ in range(NN)])
        ns = np.array([coh_score(Zc, prot, ii, rng.permutation(sv)) for _ in range(NN)]) if (sv < 0).any() and (sv > 0).any() else np.full(NN, np.nan)
        rows.append(dict(set=nm, cohort=coh, n_sites=len(sg), frac_pos=(sv > 0).mean(), coherence=obs, null_mean=np.nanmean(nm_),
                         z_matched=(obs - np.nanmean(nm_)) / np.nanstd(nm_), p_matched=(1 + (nm_ >= obs).sum()) / (1 + NN),
                         p_signshuffle=(1 + (ns >= obs).sum()) / (1 + NN) if np.isfinite(ns).any() else np.nan))
    print(coh, flush=True)
T = pd.DataFrame(rows); T.to_csv(RT + 'T8_coherence.tsv', sep='\t', index=False)
T['kind'] = T.set.str.split(':').str[0]
agg = T.groupby('set').agg(kind=('kind', 'first'), n_coh=('cohort', 'size'), n_sig=('p_matched', lambda p: (p < 0.05).sum()),
                           n_sig_sign=('p_signshuffle', lambda p: (p < 0.05).sum()), med_coh=('coherence', 'median'), med_z=('z_matched', 'median'))
agg.to_csv(RT + 'T8_coherence_summary.tsv', sep='\t')
print(agg.groupby('kind').agg(n=('n_coh', 'size'), frac_rep3=('n_sig', lambda x: (x >= 3).mean()), med_z=('med_z', 'median'), med_nsig=('n_sig', 'median')))
pd.set_option('display.width', 200); print(agg.sort_values('med_z', ascending=False).round(3).head(60).to_string())
