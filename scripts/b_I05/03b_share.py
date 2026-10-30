"""Task 1(i), noise-robust version. Per-site onset is confounded by noise: false-positive responders have
uniform 'onset' and are over-represented in the reverse direction. Population statistic instead:
reverse share at time t = n_rev(t) / (n_prim(t) + n_rev(t)), responders |log2FC|>thr at time t (sites measured at all t).
Statistic per series: slope of logit(reverse share) vs time rank (scaled 0..1, i.e. last-minus-first for 2 tp).
>0 means the reverse layer grows relative to the primary layer over time.
Also 'transient vs late-only' asymmetry: among sites responsive at exactly one time point (|x|>thr there, |x|<0.5 elsewhere)...
Controls: shuffled time labels (statistic is antisymmetric under reversal -> null centred on 0); thresholds 1 and 2.
Aggregation: series -> study -> mean; study bootstrap."""
import pickle, numpy as np, pandas as pd
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
SM = pd.read_csv(CA + 'series_meta.tsv', sep='\t', dtype={'study': str})
VV = pickle.load(open(CA + 'series_V.pkl', 'rb'))
SM['family'] = np.select([SM.prim == 'down', SM.agent == 'rituximab', SM.agent.isin(['nicotine'])], ['inhibitor', 'rituximab', 'nicotine'], 'stimulus')
SM.loc[SM.prim == 'para', 'family'] = 'paradoxical'
rng = np.random.default_rng(1)


def share_slope(V, prim, thr, perm=None):
    if perm is not None: V = V[:, perm]
    sg = 1 if prim == 'up' else -1
    npr = (sg * V > thr).sum(0) + 0.5; nrv = (sg * V < -thr).sum(0) + 0.5
    lg = np.log(nrv / npr); x = np.linspace(0, 1, V.shape[1])
    tot = npr + nrv
    if tot.min() < 5: return np.nan, lg
    return np.polyfit(x, lg, 1, w=np.sqrt(1 / (1 / npr + 1 / nrv)))[0], lg


def boot(st, nb=2000):
    st = st.dropna()
    b = [st.sample(len(st), replace=True, random_state=rng.integers(1e9)).mean() for _ in range(nb)]
    return st.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5), len(st)


rows = []
for s in SM[SM.prim.isin(['up', 'down'])].itertuples():
    V = VV[s.series][1]
    for thr in [1.0, 2.0]:
        sl, lg = share_slope(V, s.prim, thr)
        rows.append(dict(series=s.series, study=s.study, family=s.family, agent=s.agent, cell=s.cell, n_tp=s.n_tp, times=s.times, thr=thr,
                         slope=sl, logit_rev_first=lg[0], logit_rev_last=lg[-1]))
R = pd.DataFrame(rows); R.to_csv(OUT + 'T1b_revshare_per_series.tsv', sep='\t', index=False)
# null by time shuffling
NULL = {}
for it in range(300):
    for s in SM[SM.prim.isin(['up', 'down'])].itertuples():
        V = VV[s.series][1]; p = rng.permutation(V.shape[1])
        NULL.setdefault(it, []).append(dict(study=s.study, family=s.family, slope=share_slope(V, s.prim, 1.0, p)[0]))
out = []
fams = lambda d: [('inh+stim', d[d.family.isin(['inhibitor', 'stimulus'])]), ('inhibitor', d[d.family == 'inhibitor']), ('stimulus', d[d.family == 'stimulus']),
                  ('rituximab', d[d.family == 'rituximab']), ('nicotine', d[d.family == 'nicotine'])]
nulls = {}
for it, rr in NULL.items():
    rr = pd.DataFrame(rr)
    for f, d in fams(rr): nulls.setdefault(f, []).append(d.groupby('study').slope.mean().mean())
for thr, d0 in R.groupby('thr'):
    for f, d in fams(d0):
        st = d.groupby('study').slope.mean()
        m, lo, hi, ns = boot(st)
        r = dict(thr=thr, family=f, slope=m, lo=lo, hi=hi, n_studies=ns, n_series=int(d.slope.notna().sum()), frac_studies_pos=(st > 0).mean())
        if thr == 1.0:
            nn = np.array(nulls[f]); r['null_mean'] = nn.mean(); r['null_sd'] = nn.std(); r['p_perm'] = (np.sum(np.abs(nn) >= abs(m)) + 1) / (len(nn) + 1)
        out.append(r)
O = pd.DataFrame(out); O.to_csv(OUT + 'T1b_revshare_summary.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(O.round(3).to_string())
print(R[R.thr == 1].groupby(['family', 'study', 'agent']).agg(n=('series', 'size'), times=('times', 'first'), slope=('slope', 'mean'),
      first=('logit_rev_first', 'mean'), last=('logit_rev_last', 'mean')).round(2).to_string())

# ---- noise-robust: only 'confirmed' responders (same direction > 0.5 at an adjacent time point); series with >=3 time points
def share_slope_conf(V, prim, thr, perm=None):
    if perm is not None: V = V[:, perm]
    sg = 1 if prim == 'up' else -1; X = sg * V; T = X.shape[1]
    nb_p = np.zeros_like(X, bool); nb_r = np.zeros_like(X, bool)
    for t in range(T):
        nbr = [u for u in (t - 1, t + 1) if 0 <= u < T]
        nb_p[:, t] = (X[:, nbr] > 0.5).any(1); nb_r[:, t] = (X[:, nbr] < -0.5).any(1)
    npr = ((X > thr) & nb_p).sum(0) + 0.5; nrv = ((X < -thr) & nb_r).sum(0) + 0.5
    if (npr + nrv).min() < 5: return np.nan
    x = np.linspace(0, 1, T); return np.polyfit(x, np.log(nrv / npr), 1, w=np.sqrt(1 / (1 / npr + 1 / nrv)))[0]

rows = []
for s in SM[SM.prim.isin(['up', 'down']) & (SM.n_tp >= 3)].itertuples():
    V = VV[s.series][1]
    obs = share_slope_conf(V, s.prim, 1.0)
    nl = [share_slope_conf(V, s.prim, 1.0, rng.permutation(V.shape[1])) for _ in range(100)]
    rows.append(dict(series=s.series, study=s.study, family=s.family, agent=s.agent, n_tp=s.n_tp, slope_conf=obs, null_mean=np.nanmean(nl)))
Rc = pd.DataFrame(rows); Rc.to_csv(OUT + 'T1c_revshare_confirmed_per_series.tsv', sep='\t', index=False)
for f, d in fams(Rc):
    st = d.groupby('study').slope_conf.mean()
    if st.notna().sum() == 0: continue
    m, lo, hi, ns = boot(st); print('confirmed', f, round(m, 3), round(lo, 3), round(hi, 3), ns, 'null', round(d.groupby('study').null_mean.mean().mean(), 3))
print(Rc.round(2).to_string())
