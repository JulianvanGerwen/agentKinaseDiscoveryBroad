"""I10 step 3: score CPTAC tumours by each footprint (rank-based).
Per cohort & input variant (tc_padj = protein-adjusted [primary], tc = unadjusted):
  site-wise robust z across tumours -> per-sample rank -> normal scores (rank-based, removes sample-level shifts);
  score = sum w_i r_i / sum |w_i| over overlapping footprint sites (>=10 sites).
Null: NR=200 random site sets, each footprint site replaced by a random site from the same (coverage x SD) bin, keeping weights.
Outputs: cache/b_I10/cptac_scores_{variant}.parquet (long: sample, cohort, footprint, raw, nz, n_sites), cache/b_I10/cptac_null_{variant}_{cohort}.npy,
         cache/b_I10/cptac_rank_{variant}_{cohort}.parquet (rank-normal matrix, used by modules)."""
import sys, numpy as np, pandas as pd
from scipy.stats import norm, rankdata
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
NR = 200; rng = np.random.default_rng(1)
FP = pd.read_csv(RT + 'T1_footprints.tsv', sep='\t')
groups = sorted(FP.group.unique())
# also down-only / up-only variants for MAPK/PI3K/EGFR inhibitor footprints (portability test)
VAR = {}
for g in groups:
    f = FP[FP.group == g]; VAR[g] = f
    if g.startswith('KI_') or g.endswith('axis'):
        VAR[g + '|down'] = f[f.w < 0]; VAR[g + '|up'] = f[f.w > 0]
names = list(VAR)

def rank_matrix(x):
    sd = (x.sub(x.median(1), axis=0)).abs().median(1) * 1.4826
    z = x.sub(x.median(1), axis=0).div(sd.replace(0, np.nan), axis=0)
    r = z.rank(axis=0); n = z.notna().sum(0)
    return pd.DataFrame(norm.ppf((r - 0.5) / n), index=x.index, columns=x.columns).astype('float32'), sd

for variant in ['tc_padj', 'tc']:
    rows = []
    for coh in cl.COHORTS:
        x = cl.matrix(coh, f'input_{variant}')
        tum = [s for s in cl.tumour_ids(coh) if s in x.columns]
        x = x[tum]; x = x[x.notna().mean(1) >= 0.25]
        R, sd = rank_matrix(x)
        R.to_parquet(CA + f'cptac_rank_{variant}_{coh}.parquet')
        cov = x.notna().mean(1)
        b = pd.qcut(cov.rank(method='first'), 5, labels=False) * 5 + pd.qcut(sd.rank(method='first'), 5, labels=False)
        bins = {k: np.where(b.values == k)[0] for k in range(25)}
        Rv = R.values; M = np.isfinite(Rv); Rz = np.nan_to_num(Rv)
        idx = {s: i for i, s in enumerate(R.index)}
        null = np.full((len(names), NR, len(tum)), np.nan, np.float32)
        for j, nm in enumerate(names):
            f = VAR[nm]; f = f[f.site.isin(idx)]
            if len(f) < 10: continue
            ii = np.array([idx[s] for s in f.site]); w = f.w.values[:, None]
            num = (w * Rz[ii]).sum(0); den = (np.abs(w) * M[ii]).sum(0); nsite = M[ii].sum(0)
            raw = np.where(nsite >= 10, num / np.where(den > 0, den, np.nan), np.nan)
            sb = b.values[ii]
            for k in range(NR):
                jj = np.array([rng.choice(bins[q]) for q in sb])
                n2 = (w * Rz[jj]).sum(0); d2 = (np.abs(w) * M[jj]).sum(0)
                null[j, k] = n2 / np.where(d2 > 0, d2, np.nan)
            nz = (raw - np.nanmean(null[j], 0)) / np.nanstd(null[j], 0)
            rows.append(pd.DataFrame({'sample': tum, 'cohort': coh, 'footprint': nm, 'raw': raw, 'nz': nz, 'n_sites': nsite, 'n_fp_overlap': len(f)}))
        np.save(CA + f'cptac_null_{variant}_{coh}.npy', null)
        print(variant, coh, len(tum), x.shape[0], flush=True)
    S = pd.concat(rows); S.to_parquet(CA + f'cptac_scores_{variant}.parquet')
    ov = S.groupby(['footprint', 'cohort']).n_fp_overlap.first().unstack()
    ov.to_csv(RT + f'T3_cptac_overlap_{variant}.tsv', sep='\t')
pd.Series(names).to_csv(CA + 'cptac_score_names.tsv', sep='\t', index=False, header=False)
