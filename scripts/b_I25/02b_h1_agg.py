"""H1 aggregation: per study log((O+.5)/(E+.5)) averaged over studies; study-cluster bootstrap on study-level sums; LOSO; time strata."""
import numpy as np, pandas as pd
OUT = '/workspace/results/b_I25/tables/'
R = pd.read_csv(OUT + 'T1_H1_per_exp.tsv', sep='\t', dtype={'study': str})
rng = np.random.default_rng(1)
def st_tab(d):
    return d.groupby(['study', 'dir'])[['O', 'E']].sum().unstack('dir')
def lo(t, lab, w=None):
    v = np.log((t[('O', lab)] + .5) / (t[('E', lab)] + .5)); m = v.notna()
    return np.average(v[m], weights=None if w is None else w[m.values])
out = []
for (kind, sname), d in R.groupby(['kind', 'set']):
    t = st_tab(d).dropna(); ns = len(t)
    if ns < 2: continue
    rev, pri = lo(t, 'reverse'), lo(t, 'primary')
    bs = []
    for _ in range(2000):
        w = np.bincount(rng.integers(0, ns, ns), minlength=ns).astype(float)
        r_, p_ = lo(t, 'reverse', w), lo(t, 'primary', w); bs.append((r_, p_, r_ - p_))
    bs = np.array(bs); q = lambda k: np.percentile(bs[:, k], [2.5, 97.5])
    loso = [lo(t.drop(s), 'reverse') - lo(t.drop(s), 'primary') for s in t.index]
    tr = {}
    for tl, m in [('le60', d.time_min <= 60), ('gt60', d.time_min > 60)]:
        tt = st_tab(d[m]).dropna(); tr[tl] = lo(tt, 'reverse') - lo(tt, 'primary') if len(tt) >= 2 else np.nan; tr[tl + '_ns'] = len(tt)
    dr, dp = d[d.dir == 'reverse'], d[d.dir == 'primary']
    out.append(dict(kind=kind, set=sname, n_studies=ns, n_exp=d.exp.nunique(), logOE_rev=rev, rev_lo=q(0)[0], rev_hi=q(0)[1],
                    logOE_prim=pri, prim_lo=q(1)[0], prim_hi=q(1)[1], diff=rev - pri, diff_lo=q(2)[0], diff_hi=q(2)[1],
                    loso_min=min(loso), loso_max=max(loso), diff_le60=tr['le60'], ns_le60=tr['le60_ns'], diff_gt60=tr['gt60'], ns_gt60=tr['gt60_ns'],
                    rev_O=dr.O.sum(), rev_E=dr.E.sum(), prim_O=dp.O.sum(), prim_E=dp.E.sum(),
                    pooled_logOE_rev=np.log(dr.O.sum() / dr.E.sum()), pooled_logOE_prim=np.log(dp.O.sum() / dp.E.sum())))
T = pd.DataFrame(out); T.to_csv(OUT + 'T1_H1_enrichment.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); print(T.drop(columns=['rev_O', 'rev_E', 'prim_O', 'prim_E']).round(3).to_string())
# per-study detail for key sets
for s in ['CC_Mup', 'CC_union', 'CC_MKI']:
    d = R[(R.kind == 'inh') & (R.set == s)]; t = st_tab(d)
    t['n_exp'] = d.groupby('study').exp.nunique()
    t['lOE_rev'] = np.log((t[('O', 'reverse')] + .5) / (t[('E', 'reverse')] + .5)); t['lOE_prim'] = np.log((t[('O', 'primary')] + .5) / (t[('E', 'primary')] + .5))
    print(s); print(t.round(2).sort_values('n_exp', ascending=False).head(12).to_string())
