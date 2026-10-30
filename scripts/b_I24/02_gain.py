"""Per-cohort gain (robust slope) per edge; residuals; cross-cohort reproducibility; random-site null."""
from common import *
from scipy import stats
rng = np.random.default_rng(0)
COV = ['purity', 'prolif']
res, resid_all = [], []
for variant in ['tc_padj', 'tc']:
    A = load_stages(variant)
    for c in COHORTS:
        D = A[A.cohort == c]
        for e in EDGES:
            x, y = e
            d = D[[x, y] + COV].dropna()
            if len(d) < 30 or d[x].std() == 0: continue
            # raw-unit (elasticity) robust slope
            r = rlm_fit(d[y], d[[x] + COV])
            # standardized slope
            ds = d.apply(zs)
            rs = rlm_fit(ds[y], ds[[x] + COV])
            rho = stats.spearmanr(sm.OLS(stats.rankdata(d[x]), sm.add_constant(d[COV])).fit().resid,
                                  sm.OLS(stats.rankdata(d[y]), sm.add_constant(d[COV])).fit().resid)[0]
            res.append(dict(variant=variant, cohort=c, edge=ename(e), n=len(d), slope=r.params[x], se=r.bse[x],
                            std_slope=rs.params[x], std_se=rs.bse[x], partial_rho=rho))
            if variant == 'tc_padj':
                rr = pd.DataFrame({'sample_id': d.index, 'cohort': c, 'edge': ename(e), 'resid': (r.resid / r.scale).values,
                                   'x': d[x].values, 'y': d[y].values})
                resid_all.append(rr)
G = pd.DataFrame(res)
G.to_csv(f'{OUT}/tables/T2_gain_by_cohort.tsv', sep='\t', index=False)
pd.concat(resid_all).to_parquet(f'{OUT}/tables/transmission_residuals.parquet')
rows = []
for (v, e), g in G.groupby(['variant', 'edge']):
    m = meta_re(g.std_slope, g.std_se); mr = meta_re(g.slope, g.se)
    rows.append(dict(variant=v, edge=e, k=m['k'], std_slope_meta=m['est'], std_lo=m['est'] - 1.96 * m['se'], std_hi=m['est'] + 1.96 * m['se'],
                     p=m['p'], I2=m['I2'], elasticity_meta=mr['est'], el_lo=mr['est'] - 1.96 * mr['se'], el_hi=mr['est'] + 1.96 * mr['se'],
                     n_pos_sig=int(((g.std_slope / g.std_se) > 1.96).sum()), min_cohort=g.std_slope.min(), max_cohort=g.std_slope.max(),
                     median_partial_rho=g.partial_rho.median()))
M = pd.DataFrame(rows); M.to_csv(f'{OUT}/tables/T2_gain_meta.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(M.round(3).to_string())
# reproducibility of slopes across cohorts: rank correlation of per-edge slopes between cohort pairs
P = G[G.variant == 'tc_padj'].pivot(index='edge', columns='cohort', values='std_slope')
cc = P.corr(method='spearman', min_periods=6); iu = np.triu_indices_from(cc, 1)
print('median pairwise Spearman of edge-slope profiles across cohorts:', np.nanmedian(cc.values[iu]))
print(P.round(2).to_string())
