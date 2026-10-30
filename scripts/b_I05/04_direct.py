"""Task 1(ii): are early primary responders enriched for the perturbed kinase's motif / known substrates vs late ones?
Direct-target set per series: inhibitor primary targets (+ single obligate effector, e.g. BRAF->ERK1/2, MTOR->S6K1) or
stimulus_expected_up_kinases. motif-matched = max Cantley percentile over the set >= 90 (S/T sets scored on S/T sites only).
(a) per-site onset version: among primary responders (|peak|>1), MH log OR of motif-match for early-onset vs late-onset, strata=peak bin.
(b) per-time-point version: enrichment (log OR vs non-responding measured sites, same residue class) of motif-match among
    primary and among reverse responders at the FIRST vs LAST time point.
Matched random sites = the measured non-responders of the same series and residue class."""
import numpy as np, pandas as pd, pickle
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
SM = pd.read_csv(CA + 'series_meta.tsv', sep='\t', dtype={'study': str})
SM['family'] = np.select([SM.prim == 'down', SM.agent == 'rituximab', SM.agent.isin(['nicotine'])], ['inhibitor', 'rituximab', 'nicotine'], 'stimulus')
L = pd.read_parquet(CA + 'onset_long.parquet')
_S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'psp_kinases_human', 'psp_regulatory'], low_memory=False).set_index('site')
L['psp_any'] = _S.psp_kinases_human.reindex(L.site).notna().values
L['psp_reg'] = _S.psp_regulatory.reindex(L.site).fillna(False).astype(bool).values
VV = pickle.load(open(CA + 'series_V.pkl', 'rb'))
TYR = {'EGFR', 'MET', 'NTRK1', 'JAK1', 'JAK2', 'JAK3', 'TYK2', 'ZAP70', 'LCK', 'SYK', 'LYN', 'BTK'}
BINS = [1, 2, np.inf]; rng = np.random.default_rng(2)


def lor(a, b):  # a, b boolean arrays: P(a) in group vs P(b) in other group
    return np.log((a.sum() + .5) / (len(a) - a.sum() + .5)) - np.log((b.sum() + .5) / (len(b) - b.sum() + .5))


def mh(y, x, strata):
    num = den = 0.0
    for s in np.unique(strata):
        m = strata == s; n = m.sum()
        a = np.sum(y[m] & x[m]); b = np.sum(~y[m] & x[m]); c = np.sum(y[m] & ~x[m]); d = np.sum(~y[m] & ~x[m])
        num += (a + .25) * (d + .25) / n; den += (b + .25) * (c + .25) / n
    return np.log(num / den)


rows = []
for s in SM[SM.prim.isin(['up', 'down']) & SM.kin_used.notna()].itertuples():
    d = L[L.series == s.series].copy()
    ks = s.kin_used.split(';'); st_only = not (set(ks) & TYR)
    if st_only: d = d[d.res.isin(['S', 'T'])]
    d['match'] = d.tgt_pct >= 90
    sg = 1 if s.prim == 'up' else -1
    resp = d.peak.abs() > 1; prim = np.sign(d.peak) == sg
    P = d[resp & prim]; Rv = d[resp & ~prim]; bg = d[~resp]
    early = (P["on1.0"] == 0).values; strata = np.digitize(P.peak.abs(), BINS)
    r = dict(series=s.series, study=s.study, family=s.family, agent=s.agent, kin=s.kin_used, n_tp=s.n_tp, n_prim=len(P), n_rev=len(Rv),
             bg_match=d.match.mean(), prim_match=P.match.mean(), rev_match=Rv.match.mean(),
             lor_prim_vs_bg=lor(P.match.values, bg.match.values), lor_rev_vs_bg=lor(Rv.match.values, bg.match.values),
             lor_psp_prim_vs_bg=lor(P.tgt_psp.values, bg.tgt_psp.values), lor_psp_rev_vs_bg=lor(Rv.tgt_psp.values, bg.tgt_psp.values))
    if early.sum() >= 5 and (~early).sum() >= 5:
        r['mh_match_early_vs_late'] = mh(P.match.values, early, strata)
        r['mh_psp_early_vs_late'] = mh(P.tgt_psp.values, early, strata)
        nt = ~P.tgt_psp.values
        r['mh_pspOTHER_early_vs_late'] = mh((P.psp_any.values & nt)[nt], early[nt], strata[nt])
        r['mh_pspreg_early_vs_late'] = mh(P.psp_reg.values, early, strata)
        Rv_e = (Rv['on1.0'] == 0).values
        if Rv_e.sum() >= 5 and (~Rv_e).sum() >= 5: r['mh_psp_REV_early_vs_late'] = mh(Rv.tgt_psp.values, Rv_e, np.digitize(Rv.peak.abs(), BINS)); r['mh_match_REV_early_vs_late'] = mh(Rv.match.values, Rv_e, np.digitize(Rv.peak.abs(), BINS))
        r['n_early'] = int(early.sum()); r['n_late'] = int((~early).sum())
        r['match_early'] = P.match.values[early].mean(); r['match_late'] = P.match.values[~early].mean()
    # per time point (b)
    sites, V, times = VV[s.series]
    idx = pd.Index(sites).get_indexer(d.site); Vd = V[idx] * sg; mt = d.match.values
    for lab, t in [('first', 0), ('last', V.shape[1] - 1)]:
        pr = Vd[:, t] > 1; rv = Vd[:, t] < -1; nr = np.abs(Vd).max(1) <= 1
        r[f'tp_{lab}_lor_prim'] = lor(mt[pr], mt[nr]) if pr.sum() >= 5 else np.nan
        r[f'tp_{lab}_lor_rev'] = lor(mt[rv], mt[nr]) if rv.sum() >= 5 else np.nan
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(OUT + 'T2_directness_per_series.tsv', sep='\t', index=False)


def boot(st, nb=2000):
    st = st.dropna()
    if len(st) == 0: return [np.nan] * 4
    b = [st.sample(len(st), replace=True, random_state=rng.integers(1e9)).mean() for _ in range(nb)]
    return st.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5), len(st)


R['tp_delta_prim'] = R.tp_first_lor_prim - R.tp_last_lor_prim
R['tp_delta_rev'] = R.tp_first_lor_rev - R.tp_last_lor_rev
R['prim_minus_rev'] = R.lor_prim_vs_bg - R.lor_rev_vs_bg
out = []
for f, d in [('inh+stim', R[R.family.isin(['inhibitor', 'stimulus'])]), ('inhibitor', R[R.family == 'inhibitor']), ('stimulus', R[R.family == 'stimulus'])]:
    for col in ['lor_prim_vs_bg', 'lor_rev_vs_bg', 'prim_minus_rev', 'lor_psp_prim_vs_bg', 'lor_psp_rev_vs_bg', 'mh_match_early_vs_late', 'mh_psp_early_vs_late', 'mh_pspOTHER_early_vs_late', 'mh_pspreg_early_vs_late', 'mh_psp_REV_early_vs_late', 'mh_match_REV_early_vs_late',
                'tp_first_lor_prim', 'tp_last_lor_prim', 'tp_delta_prim', 'tp_first_lor_rev', 'tp_last_lor_rev', 'tp_delta_rev']:
        m, lo, hi, ns = boot(d.groupby('study')[col].mean())
        out.append(dict(family=f, metric=col, est=m, lo=lo, hi=hi, n_studies=ns))
O = pd.DataFrame(out); O.to_csv(OUT + 'T2_directness_summary.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(O.round(3).to_string())
print(R.groupby(['family', 'study', 'agent'])[['n_prim', 'n_rev', 'bg_match', 'prim_match', 'rev_match', 'lor_prim_vs_bg', 'lor_rev_vs_bg', 'mh_match_early_vs_late', 'tp_delta_prim']].mean().round(2).to_string())
