"""Shared helpers for the I19 audit (independent re-implementation of the replication metric)."""
import numpy as np, pandas as pd
CA = "/workspace/results/audit_I19/cache/"
BINS = [(0.5, 1), (1, 2), (2, 99)]
V600 = {'COLO205', 'WIDR', 'A375', 'WM239A', 'HT29'}


def load_exp(name='inh'):
    M = pd.read_csv(CA + f'{name}_meta.tsv', sep='\t', dtype={'study_pmid': str}).fillna({'binders': ''})
    X = pd.read_parquet(CA + f'{name}_exp.parquet')
    E = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), agent=('agent', 'first'), cls=('cls', 'first'),
                             binders=('binders', 'first'), time_min=('time_min', 'median'), tmin=('time_min', 'min'), tmax=('time_min', 'max'),
                             dose_nM=('dose_nM', 'median'), n_cond=('conditionID', 'size')).reset_index()
    E = E.set_index('exp').loc[X.columns].reset_index()
    return X, E


def mask_binders(X, E):
    S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'gene'], low_memory=False).set_index('site')
    gene = S.gene.reindex(X.index).fillna('').values
    A = X.values.astype(np.float32).copy()
    for j, b in enumerate(E.binders):
        bs = set(x.strip() for x in str(b).replace('/', ';').split(';') if x.strip())
        if bs: A[np.isin(gene, list(bs)), j] = np.nan
    return A


def center(A, E, min_exp=5, min_cls=3):
    """subtract per-site median across experiments of a study x cell group (>=5 exps, >=3 classes); others unchanged."""
    A = A.copy(); g = E.study + '|' + E.cell.astype(str)
    for k, d in E.groupby(g):
        if len(d) >= min_exp and d.cls.nunique() >= min_cls:
            ix = d.index.values
            with np.errstate(all='ignore'):
                A[:, ix] = A[:, ix] - np.nanmedian(A[:, ix], 1, keepdims=True)
    return A


def pair_rows(A, E, groups, thr=0.3, bins=BINS, min_disc=5, site_mask=None, pair_filter=None, cross_study=True, null_pool='other_class_same_study', min_both=200):
    """Experiment-pair replication. For ordered pairs (i,j) within a group (different studies if cross_study):
    sites measured in both (and in site_mask); discovery: dir*x_i in bin; replication: dir*x_j > thr.
    base = P(dir*x_j > thr) over all co-measured sites (per-direction, per-experiment base rate; hypergeometric expectation);
    snull = mean over k in null pool (other-class exps of study_j) of P(dir*x_k > thr | discovery sites measured in k).
    lor = log odds ratio (disc vs non-disc co-measured) x (rep vs not)."""
    rows = []
    studies = E.study.values; cls = E.cls.values
    finite = ~np.isnan(A)
    for g, idx in groups.items():
        ids = np.where(idx)[0]
        for i in ids:
            for j in ids:
                if i == j: continue
                if cross_study and studies[i] == studies[j]: continue
                if pair_filter is not None and not pair_filter(i, j): continue
                both = finite[:, i] & finite[:, j]
                if site_mask is not None: both &= site_mask
                if both.sum() < min_both: continue
                ai, aj = A[both, i], A[both, j]
                pool = []
                if null_pool == 'other_class_same_study':
                    pool = np.where((studies == studies[j]) & (cls != cls[j]) & ~idx)[0]
                for sgn, lab in [(1, 'up'), (-1, 'down')]:
                    x, y = sgn * ai, sgn * aj
                    rep_all = y > thr; base = rep_all.mean()
                    for lo, hi in bins:
                        d = (x > lo) & (x <= hi)
                        nd = int(d.sum())
                        if nd < min_disc: continue
                        r = rep_all[d].mean(); rn = rep_all[~d].mean()
                        a_, b_, c_, d_ = rep_all[d].sum() + .5, nd - rep_all[d].sum() + .5, rep_all[~d].sum() + .5, (~d).sum() - rep_all[~d].sum() + .5
                        row = dict(group=g, i=i, j=j, si=studies[i], sj=studies[j], dir=lab, bin=f'{lo}-{hi}', n_disc=nd, n_both=int(both.sum()),
                                   rep=r, base=base, lor=np.log(a_ * d_ / (b_ * c_)), meanB=float(y[d].mean()), opp=float(np.mean(y[d] < -thr)), opp_base=float(np.mean(y < -thr)))
                        if len(pool):
                            sites_idx = np.where(both)[0][d]
                            vals = []
                            for k in pool:
                                v = sgn * A[sites_idx, k]; v = v[~np.isnan(v)]
                                if len(v) >= 5: vals.append(np.mean(v > thr))
                            row['snull'] = np.mean(vals) if vals else np.nan
                        rows.append(row)
    return pd.DataFrame(rows)


def summarise(R, by=('dir', 'bin'), weight='study_pair'):
    """weight='sites' (pool discovery sites, like I19) or 'study_pair' (mean over pairs within study pair, then over study pairs)."""
    by = list(by); out = []
    for k, d in R.groupby(by):
        r = dict(zip(by, k if isinstance(k, tuple) else (k,)))
        r['pairs'] = len(d); r['study_pairs'] = d.groupby(['si', 'sj']).ngroups; r['n_disc'] = int(d.n_disc.sum())
        cols = [c for c in ['rep', 'base', 'snull', 'lor', 'meanB', 'opp', 'opp_base'] if c in d]
        if weight == 'sites':
            for c in cols:
                m = d[c].notna(); r[c] = np.average(d.loc[m, c], weights=d.loc[m, 'n_disc']) if m.any() else np.nan
        else:
            sp = d.groupby(['si', 'sj'])[cols].mean()
            for c in cols: r[c] = sp[c].mean()
            # bootstrap over studies (cluster): resample study pairs
        out.append(r)
    return pd.DataFrame(out)


def boot_diff(R, metric='lor', bin_=None, nboot=500, seed=0):
    """cluster bootstrap over studies of (down - up) difference in study-pair-mean metric."""
    rng = np.random.default_rng(seed)
    d = R if bin_ is None else R[R.bin == bin_]
    sp = d.groupby(['si', 'sj', 'dir'])[metric].mean().unstack('dir').dropna()
    st = np.unique(np.r_[sp.index.get_level_values(0), sp.index.get_level_values(1)])
    obs = (sp['down'] - sp['up']).mean(); bs = []
    si = sp.index.get_level_values(0).values; sj = sp.index.get_level_values(1).values; diff = (sp['down'] - sp['up']).values
    for _ in range(nboot):
        pick = rng.choice(st, len(st), replace=True); w = pd.Series(pick).value_counts()
        ww = w.reindex(si).fillna(0).values * w.reindex(sj).fillna(0).values
        if ww.sum() > 0: bs.append(np.average(diff, weights=ww))
    return obs, np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(st), len(sp)


def load_cond(name='inh'):
    M = pd.read_csv(CA + f'{name}_meta.tsv', sep='\t', dtype={'study_pmid': str}).fillna({'binders': ''})
    X = pd.read_parquet(CA + f'{name}_cond.parquet')[M.conditionID]
    E = M.rename(columns={'study_pmid': 'study', 'cell_line_norm': 'cell'}).reset_index(drop=True)
    return X, E


def add_opp(R, A, E, thr=0.3):
    """not used"""
    return R
