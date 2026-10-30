"""Task 1 + controls.
(i) Are primary-direction responders earlier than reverse-direction responders at matched peak magnitude?
    per series: Mantel-Haenszel log OR of 'early' (onset at first time point) for primary vs reverse, strata = |peak| bin.
    Also mean normalised onset (onset index / (n_tp-1)) difference, peak-bin matched (multi-time-point series).
(ii) Directness: among primary responders, target-kinase motif match (Cantley pct >= 90) and PSP known target substrates:
    MH log OR (direct | early vs late), strata = peak bin. Also each class vs residue-matched random measured sites.
Aggregation: series -> study (mean) -> across studies; 95% CI by study bootstrap (2000).
Controls: shuffled time labels within series (same permutation for all sites; 200 perms); thresholds 0.5 / 1.0 / half-max."""
import pickle, numpy as np, pandas as pd
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
SM = pd.read_csv(CA + 'series_meta.tsv', sep='\t', dtype={'study': str})
L = pd.read_parquet(CA + 'onset_long.parquet')
VV = pickle.load(open(CA + 'series_V.pkl', 'rb'))
SM['family'] = np.select([SM.prim == 'down', SM.agent == 'rituximab', SM.agent.isin(['nicotine'])], ['inhibitor', 'rituximab', 'nicotine'], 'stimulus')
SM.loc[SM.prim == 'para', 'family'] = 'paradoxical'
BINS = [0.5, 1, 2, np.inf]
rng = np.random.default_rng(0)


def mh_lor(y, x, strata):
    """Mantel-Haenszel log OR of y (bool) for x (bool) vs ~x, stratified."""
    num = den = 0.0
    for s in np.unique(strata):
        m = strata == s; n = m.sum()
        if n < 4: continue
        a = np.sum(y[m] & x[m]); b = np.sum(~y[m] & x[m]); c = np.sum(y[m] & ~x[m]); d = np.sum(~y[m] & ~x[m])
        num += (a + .25) * (d + .25) / n; den += (b + .25) * (c + .25) / n
    return np.log(num / den) if den > 0 and num > 0 else np.nan


def onset_from(V, thr, perm=None):
    if perm is not None: V = V[:, perm]
    i = np.abs(V).argmax(1); pk = V[np.arange(len(V)), i]; dv = V * np.sign(pk)[:, None]
    if thr == 'half': hit = dv >= 0.5 * np.abs(pk)[:, None]
    else: hit = dv > thr
    return pk, np.where(hit.any(1), hit.argmax(1), -1)


def series_stats(V, prim, thr, perm=None, min_n=5):
    pk, on = onset_from(V, thr, perm)
    t = 1.0 if thr == 'half' else thr
    resp = np.abs(pk) > max(t, 0.5) if thr != 'half' else np.abs(pk) > 1.0
    if resp.sum() < 10: return None
    p = (np.sign(pk) == (1 if prim == 'up' else -1))[resp]
    if p.sum() < min_n or (~p).sum() < min_n: return None
    early = (on[resp] == 0)
    strata = np.digitize(np.abs(pk[resp]), BINS)
    onorm = on[resp] / (V.shape[1] - 1)
    # peak-matched normalised onset difference: mean within stratum, weighted by min count
    dd, ww = [], []
    for s in np.unique(strata):
        m = strata == s
        if (p & m).sum() >= 3 and (~p & m).sum() >= 3:
            dd.append(onorm[m & ~p].mean() - onorm[m & p].mean()); ww.append(min((p & m).sum(), (~p & m).sum()))
    return dict(lor_early=mh_lor(early, p, strata), d_onset=np.average(dd, weights=ww) if dd else np.nan,
                n_prim=int(p.sum()), n_rev=int((~p).sum()), early_prim=early[p].mean(), early_rev=early[~p].mean(), frac_prim=p.mean())


def boot(df, col, nb=2000):
    """series -> study mean -> mean over studies; CI from study bootstrap."""
    st = df.groupby('study')[col].mean().dropna()
    if len(st) == 0: return np.nan, np.nan, np.nan, 0
    b = [st.sample(len(st), replace=True, random_state=rng.integers(1e9)).mean() for _ in range(nb)]
    return st.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5), len(st)


rows = []
for thr in [0.5, 1.0, 'half']:
    for s in SM.itertuples():
        if s.prim not in ('up', 'down'): continue
        sites, V, times = VV[s.series]
        r = series_stats(V, s.prim, thr)
        if r is None: continue
        rows.append(dict(series=s.series, study=s.study, family=s.family, rcls=s.rcls, agent=s.agent, cell=s.cell, n_tp=s.n_tp, t0=times[0], thr=str(thr), **r))
R = pd.DataFrame(rows)
R.to_csv(OUT + 'T1_timing_per_series.tsv', sep='\t', index=False)
# shuffled-time null (aggregate statistic over all eligible series, families pooled & per family)
null = []
for it in range(200):
    rr = []
    for s in SM.itertuples():
        if s.prim not in ('up', 'down'): continue
        V = VV[s.series][1]; perm = rng.permutation(V.shape[1])
        r = series_stats(V, s.prim, 1.0, perm)
        if r: rr.append(dict(study=s.study, family=s.family, lor_early=r['lor_early'], d_onset=r['d_onset']))
    rr = pd.DataFrame(rr)
    for fam, d in [('all_excl_ritux_nic', rr[~rr.family.isin(['rituximab', 'nicotine'])])] + list(rr.groupby('family')):
        null.append(dict(it=it, family=fam, lor_early=d.groupby('study').lor_early.mean().mean(), d_onset=d.groupby('study').d_onset.mean().mean()))
null = pd.DataFrame(null)
out = []
for thr, d0 in R.groupby('thr'):
    for fam, d in [('all_excl_ritux_nic', d0[~d0.family.isin(['rituximab', 'nicotine'])])] + list(d0.groupby('family')):
        for col in ['lor_early', 'd_onset']:
            m, lo, hi, ns = boot(d, col)
            r = dict(thr=thr, family=fam, metric=col, est=m, lo=lo, hi=hi, n_studies=ns, n_series=int(d[col].notna().sum()),
                     early_prim=d.groupby('study').early_prim.mean().mean(), early_rev=d.groupby('study').early_rev.mean().mean())
            if thr == '1.0':
                nn = null[(null.family == fam)][col]
                r['null_mean'] = nn.mean(); r['null_sd'] = nn.std(); r['p_perm'] = (np.sum(np.abs(nn - nn.mean()) >= abs(m - nn.mean())) + 1) / (len(nn) + 1)
            out.append(r)
O = pd.DataFrame(out)
O.to_csv(OUT + 'T1_timing_summary.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(O.round(3).to_string())
print(R[R.thr == '1.0'].groupby(['family', 'study', 'agent']).agg(n=('series', 'size'), lor=('lor_early', 'mean'), d_on=('d_onset', 'mean'),
      ep=('early_prim', 'mean'), er=('early_rev', 'mean'), fp=('frac_prim', 'mean'), nprim=('n_prim', 'sum'), nrev=('n_rev', 'sum')).round(2).to_string())
