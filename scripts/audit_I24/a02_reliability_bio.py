"""(a) strong-instrument metas; (b) tc protein-abundance dissection; (c) split-half reliabilities, within-adaptor IV
(reliability-robust comparator) and simulation null under full transitivity with observed reliabilities;
(d) what explains RTK-independent adaptor pY (SFK substrates, immune/stroma, purity, bgPY, PTPRC)."""
from common_a import *
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
A = pd.read_parquet(f'{OUT}/tables/stages_audit.parquet')
R1 = pd.read_csv(f'{OUT}/tables/A1_decomp_by_cohort.tsv', sep='\t')
out = []
# (a) strong instruments only
for a in ['main', 'no_amp', 'EGFR_only', 'EGFR_only_no_amp', 'adj_bgPY', 'tc_adj_protein']:
    for v in ['tc', 'tc_padj']:
        g = R1[(R1.analysis == a) & (R1.variant == v) & (R1.F >= 10)]
        if len(g) < 2: continue
        m1, m2, m3 = meta_re(g.iv, g.iv_se), meta_re(g.other, g.other_se), meta_re(g['diff'], g.diff_se)
        out.append(dict(analysis=a + '_F>=10', variant=v, k=len(g), cohorts=','.join(g.cohort), iv=m1['est'], iv_lo=m1['lo'], iv_hi=m1['hi'],
                        other=m2['est'], other_lo=m2['lo'], other_hi=m2['hi'], diff=m3['est'], diff_p=m3['p']))
S1 = pd.DataFrame(out); print(S1.round(3).to_string())

# (b) tc dissection: which protein covariate removes the gap
rows = []
for c in COH:
    D = A[(A.variant == 'tc') & (A.cohort == c)]
    for name, cov in [('tc_adj_RTKprot', ('purity', 'prolif', 'RTK_prot')), ('tc_adj_ADAPT_ERKprot', ('purity', 'prolif', 'ADAPT_prot', 'ERK_prot')),
                      ('tc_adj_ERKprot', ('purity', 'prolif', 'ERK_prot')), ('tc_adj_ADAPTprot', ('purity', 'prolif', 'ADAPT_prot'))]:
        r = decomp(D, cov=cov)
        if r: r.update(cohort=c, analysis=name, variant='tc'); rows.append(r)
R2 = pd.DataFrame(rows); R2.to_csv(f'{OUT}/tables/A2_protein_dissection_by_cohort.tsv', sep='\t', index=False)
M2 = meta_table(R2, ['analysis', 'variant']); M2.to_csv(f'{OUT}/tables/A2_protein_dissection_meta.tsv', sep='\t', index=False)
print(M2[['analysis', 'k', 'fs', 'rf', 'iv', 'iv_lo', 'iv_hi', 'other', 'other_lo', 'other_hi', 'diff', 'diff_p']].round(3).to_string())
# correlation of abundance with phospho stages (tc) - how much of tc stage is abundance
ab = []
for c in COH:
    D = A[(A.variant == 'tc') & (A.cohort == c)]
    for st, p in [('RTK', 'RTK_prot'), ('EGFRpY', 'EGFR_prot'), ('ADAPT', 'ADAPT_prot'), ('ERK', 'ERK_prot'), ('ADAPT_prot', 'ERK_prot'), ('RTK_prot', 'ERK_prot'), ('RTK_prot', 'ADAPT_prot')]:
        d = D[[st, p]].dropna()
        if len(d) > 30: ab.append(dict(cohort=c, pair=f'{st}~{p}', r=stats.spearmanr(d[st], d[p])[0], n=len(d)))
AB = pd.DataFrame(ab).pivot(index='cohort', columns='pair', values='r'); AB.to_csv(f'{OUT}/tables/A2_abundance_corr.tsv', sep='\t'); print(AB.round(2).to_string())

# (c) reliability
def resid(v, C):
    C1 = np.c_[np.ones(len(v)), C]; return v - C1 @ np.linalg.lstsq(C1, v, rcond=None)[0]
def sb(a, b):
    r = np.corrcoef(a, b)[0, 1]; return 2 * r / (1 + r) if r > 0 else np.nan, r
rel, ivs = [], []
for v in ['tc', 'tc_padj']:
    for c in COH:
        X = cl.matrix(c, f'input_{v}'); D = A[(A.variant == v) & (A.cohort == c)]; X = X[D.index]
        H = {}
        rs = [s for s in SITES['RTK'] if s in X.index]
        if len(rs) >= 2: H['RTK'] = (X.loc[rs[0::2]].mean(0), X.loc[rs[1::2]].mean(0))
        H['ADAPT'] = (D.ADAPT_SHC1, D.ADAPT_PG); H['ERK'] = (D.ERK3, D.ERK2)
        for st, (h1, h2) in H.items():
            d = pd.DataFrame({'a': h1, 'b': h2, 'p': D.purity, 'q': D.prolif}).dropna()
            if len(d) < 30: continue
            C = d[['p', 'q']].values; rr, r = sb(resid(d.a.values, C), resid(d.b.values, C))
            rel.append(dict(variant=v, cohort=c, stage=st, n=len(d), half_r=r, rel_SB=rr))
        # within-adaptor IV: ERK on ADAPT_PG instrumented by ADAPT_SHC1, and vice versa; plus OLS (adj purity/prolif)
        for m, z in [('ADAPT_PG', 'ADAPT_SHC1'), ('ADAPT_SHC1', 'ADAPT_PG')]:
            d = D[[m, z, 'ERK', 'purity', 'prolif']].dropna()
            if len(d) < 30: continue
            Z = d.apply(zs); C = np.c_[np.ones(len(d)), Z[['purity', 'prolif']].values]
            def f(i):
                fs = np.linalg.lstsq(np.c_[Z[z].values[i], C[i]], Z[m].values[i], rcond=None)[0][0]
                rf = np.linalg.lstsq(np.c_[Z[z].values[i], C[i]], Z.ERK.values[i], rcond=None)[0][0]
                ols = np.linalg.lstsq(np.c_[Z[m].values[i], C[i]], Z.ERK.values[i], rcond=None)[0][0]
                return np.array([fs, rf / fs, ols])
            e = f(np.arange(len(d))); rng = np.random.default_rng(1); B = np.array([f(rng.integers(0, len(d), len(d))) for _ in range(1000)])
            se = (np.percentile(B, 75, 0) - np.percentile(B, 25, 0)) / 1.349
            ivs.append(dict(variant=v, cohort=c, mediator=m, instrument=z, n=len(d), fs=e[0], iv=e[1], iv_se=se[1], ols=e[2], ols_se=se[2]))
REL = pd.DataFrame(rel); REL.to_csv(f'{OUT}/tables/A2_split_half_reliability.tsv', sep='\t', index=False)
print(REL.pivot_table(index=['variant', 'cohort'], columns='stage', values='rel_SB').round(2).to_string())
IV = pd.DataFrame(ivs); IV.to_csv(f'{OUT}/tables/A2_within_adaptor_iv.tsv', sep='\t', index=False)
for (v, m), g in IV.groupby(['variant', 'mediator']):
    rtk = g.cohort.isin(R1[R1.analysis == 'main'].cohort.unique())
    for lab, gg in [('all10', g), ('7RTKcohorts', g[rtk])]:
        a, b = meta_re(gg.iv, gg.iv_se), meta_re(gg.ols, gg.ols_se)
        print(f'within-adaptor IV {v} {m} [{lab}] k={a["k"]}: IV {a["est"]:.3f} [{a["lo"]:.3f},{a["hi"]:.3f}]  OLS {b["est"]:.3f}')

# simulation null: full transitivity (ERK depends on adaptor only), observed n, FS, total gain and split-half reliabilities
sim = []
rng = np.random.default_rng(7)
main = R1[R1.analysis == 'main']
for v in ['tc', 'tc_padj']:
    for _, r in main[main.variant == v].iterrows():
        rl = REL[(REL.variant == v) & (REL.cohort == r.cohort)].set_index('stage').rel_SB
        relR, relA, relE = [float(np.clip(rl.get(s, 0.7), 0.2, 0.98)) for s in ['RTK', 'ADAPT', 'ERK']]
        n = int(r.n)
        for rep in range(300):
            # latent: rho = true RTK, a = lam*rho + u, e = beta*a + w ; choose lam, beta so observed fs, total roughly match
            lam = min(0.95, r.fs / np.sqrt(relR * relA)); beta = min(0.95, max(0.05, r.total / np.sqrt(relA * relE)))
            rho = rng.standard_normal(n); a = lam * rho + np.sqrt(1 - lam ** 2) * rng.standard_normal(n)
            e = beta * a + np.sqrt(1 - beta ** 2) * rng.standard_normal(n)
            def obs(x, rel): return x + np.sqrt((1 - rel) / rel) * rng.standard_normal(n)
            D = pd.DataFrame({'RTK': obs(rho, relR), 'ADAPT': obs(a, relA), 'ERK': obs(e, relE),
                              'purity': rng.standard_normal(n), 'prolif': rng.standard_normal(n)})
            q = decomp(D, nboot=1)
            sim.append(dict(variant=v, cohort=r.cohort, rep=rep, relR=relR, relA=relA, relE=relE, iv=q['iv'], other=q['other'], diff=q['diff'], F=q['F']))
SIM = pd.DataFrame(sim); SIM.to_parquet(f'{OUT}/tables/A2_sim_transitive_null.parquet')
sm_ = SIM.groupby(['variant', 'cohort']).agg(relR=('relR', 'first'), relA=('relA', 'first'), relE=('relE', 'first'), iv_med=('iv', 'median'),
                                             other_med=('other', 'median'), diff_med=('diff', 'median'), frac_diff_lt_obs=('diff', lambda x: np.nan))
for (v, c), g in SIM.groupby(['variant', 'cohort']):
    od = main[(main.variant == v) & (main.cohort == c)]['diff'].iloc[0]; sm_.loc[(v, c), 'obs_diff'] = od; sm_.loc[(v, c), 'frac_diff_lt_obs'] = np.mean(g['diff'] <= od)
# pooled null: average of 7 per-cohort diffs
for v in ['tc', 'tc_padj']:
    P = SIM[SIM.variant == v].pivot(index='rep', columns='cohort', values='diff').mean(1); od = main[main.variant == v]['diff'].mean()
    print(f'SIM {v}: null mean-of-cohort diff median {P.median():.3f} [2.5%,97.5%] [{P.quantile(.025):.3f},{P.quantile(.975):.3f}]; observed {od:.3f}; P(null<=obs)={np.mean(P <= od):.4f}')
sm_.to_csv(f'{OUT}/tables/A2_sim_summary.tsv', sep='\t'); print(sm_.round(3).to_string())

# (d) what explains RTK-independent adaptor pY
rows = []
for v in ['tc', 'tc_padj']:
    for c in COH:
        D = A[(A.variant == v) & (A.cohort == c)].copy()
        base = ['purity', 'prolif'] + (['RTK'] if D.RTK.notna().mean() > 0.4 else [])
        d = D[['ADAPT', 'ERK'] + base].dropna()
        Z = d.apply(zs); Cm = np.c_[np.ones(len(d)), Z[base].values]
        ra = Z.ADAPT.values - Cm @ np.linalg.lstsq(Cm, Z.ADAPT.values, rcond=None)[0]
        rr = pd.Series(ra, index=d.index)
        for k in ['SFK', 'Immune', 'Stromal', 'purity', 'bgPY', 'PTPRC_prot', 'ADAPT_prot', 'ERK_prot']:
            x = D[k].reindex(d.index); ok = x.notna()
            if ok.sum() < 30 or x[ok].std() == 0: continue
            rows.append(dict(variant=v, cohort=c, var=k, rho=stats.spearmanr(rr[ok], x[ok])[0], n=int(ok.sum()), rtk_in_model='RTK' in base))
        # variance of ADAPT residual explained by SFK+immune+stroma+bgPY
        ks = [k for k in ['SFK', 'Immune', 'Stromal', 'bgPY'] if D[k].notna().mean() > 0.9]
        dd = D[ks].reindex(d.index).dropna(); y = rr.reindex(dd.index).values; Xk = np.c_[np.ones(len(dd)), dd.apply(zs).values]
        b = np.linalg.lstsq(Xk, y, rcond=None)[0]; r2 = 1 - np.var(y - Xk @ b) / np.var(y)
        rows.append(dict(variant=v, cohort=c, var='R2_' + '+'.join(ks), rho=r2, n=len(dd), rtk_in_model='RTK' in base))
B = pd.DataFrame(rows); B.to_csv(f'{OUT}/tables/A2_adaptor_residual_correlates.tsv', sep='\t', index=False)
print(B[~B['var'].str.startswith('R2')].pivot_table(index=['variant', 'cohort'], columns='var', values='rho').round(2).to_string())
print(B[B['var'].str.startswith('R2')][['variant', 'cohort', 'var', 'rho']].round(3).to_string())
