"""(a) Non-transitivity: is RTK-driven adaptor pY transmitted to ERK? ERK ~ ADAPT_rtk(fitted) + ADAPT_other(resid);
(b) GBM EGFR-amp / NF1: slopes and offsets in groups using unadjusted (tc) EGFR pY;
(c) 17q11/NF1 CNV x RTK>ADAPT: robust model, PCs, leave-one-cohort-out."""
from common import *
from scipy import stats
COV = ['purity', 'prolif']; pd.set_option('display.width', 250)
rows = []
for variant in ['tc_padj', 'tc']:
    A = load_stages(variant)
    for c in COHORTS:
        D = A[A.cohort == c][['RTK', 'ADAPT', 'ERK', 'ERKSUB', 'AKTSUB'] + COV].dropna()
        if len(D) < 30: continue
        Z = D.apply(zs); C = sm.add_constant(Z[COV])
        f = sm.OLS(Z.ADAPT, sm.add_constant(Z[['RTK'] + COV])).fit()
        a_rtk = f.params['RTK'] * Z.RTK; a_oth = f.resid
        for out in ['ERK', 'ERKSUB', 'AKTSUB']:
            m = sm.RLM(Z[out], sm.add_constant(pd.DataFrame({'a_rtk': a_rtk, 'a_oth': a_oth, 'purity': Z.purity, 'prolif': Z.prolif})),
                       M=sm.robust.norms.HuberT()).fit()
            rows.append(dict(variant=variant, cohort=c, out=out, n=len(D), b_rtk_part=m.params['a_rtk'], se_rtk=m.bse['a_rtk'],
                             b_other_part=m.params['a_oth'], se_other=m.bse['a_oth']))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/tables/T7_adaptor_decomposition.tsv', sep='\t', index=False)
o = []
for (v, out), d in R.groupby(['variant', 'out']):
    a = meta_re(d.b_rtk_part, d.se_rtk); b = meta_re(d.b_other_part, d.se_other)
    d2 = (d.b_rtk_part - d.b_other_part); s2 = np.sqrt(d.se_rtk**2 + d.se_other**2); dd = meta_re(d2, s2)
    o.append(dict(variant=v, out=out, k=a['k'], rtk_part=a['est'], rtk_lo=a['est'] - 1.96 * a['se'], rtk_hi=a['est'] + 1.96 * a['se'],
                  other_part=b['est'], oth_lo=b['est'] - 1.96 * b['se'], oth_hi=b['est'] + 1.96 * b['se'], diff=dd['est'], diff_p=dd['p']))
O = pd.DataFrame(o); O.to_csv(f'{OUT}/tables/T7_adaptor_decomposition_meta.tsv', sep='\t', index=False); print(O.round(4).to_string())

print('NOTE: parts in ADAPT SD units; rtk_part = 2SLS (RTK-instrumented) gain, other_part = gain of RTK-independent adaptor pY')
# (b) GBM
for variant in ['tc', 'tc_padj']:
    A = load_stages(variant); D = A[A.cohort == 'GBM'].copy()
    gi = cl.matrix('GBM', 'cnv_gistic'); mb = cl.matrix('GBM', 'mut_binary')
    D['EGFRamp'] = (gi.loc['EGFR'].reindex(D.case_id).values == 2).astype(int)
    D['NF1mut'] = mb.loc['NF1'].reindex(D.case_id).fillna(0).values.astype(int)
    D['grp'] = np.where(D.EGFRamp == 1, 'EGFRamp', np.where(D.NF1mut == 1, 'NF1mut', 'other'))
    print(variant, D.grp.value_counts().to_dict())
    for e in [('EGFRpY', 'ADAPT'), ('ADAPT', 'ERK'), ('EGFRpY', 'ERK'), ('ERK', 'ERKSUB'), ('EGFRpY', 'ERKSUB')]:
        x, y = e; d = D[[x, y, 'grp', 'purity', 'prolif']].dropna()
        line = f'  {x}>{y}: '
        for g in ['EGFRamp', 'NF1mut', 'other']:
            q = d[d.grp == g]; mx = q[x].mean(); my = q[y].mean()
            sl = sm.RLM(q[y], sm.add_constant(q[[x]]), M=sm.robust.norms.HuberT()).fit().params[x] if len(q) >= 8 else np.nan
            line += f'{g}(n={len(q)}) mean_x={mx:.2f} mean_y={my:.2f} slope={sl:.2f} | '
        # interaction amp vs rest, standardized, with offset
        dd = d.copy(); dd['amp'] = (dd.grp == 'EGFRamp').astype(float); xz = zs(dd[x])
        M = sm.add_constant(np.column_stack([xz, dd.amp, xz * dd.amp, zs(dd.purity), zs(dd.prolif)]))
        r = sm.RLM(zs(dd[y]), M, M=sm.robust.norms.HuberT()).fit()
        line += f'amp offset {r.params[2]:.2f} (p {2*stats.norm.sf(abs(r.params[2]/r.bse[2])):.2g}), amp x slope {r.params[3]:.2f} (p {2*stats.norm.sf(abs(r.params[3]/r.bse[3])):.2g})'
        print(line)

# (c) 17q11 / NF1 CNV robustness on RTK>ADAPT
A = load_stages('tc_padj'); PCS = pd.read_parquet(f'{OUT}/tables/phospho_pcs.parquet'); rows = []
for g in ['NF1', 'BLMH']:
    for c in COHORTS:
        D = A[A.cohort == c][['RTK', 'ADAPT'] + COV + ['case_id']].dropna()
        if len(D) < 30: continue
        l = cl.matrix(c, 'cnv_log2').loc[g].reindex(D.case_id).values.astype(float); D = D[np.isfinite(l)]; l = l[np.isfinite(l)]
        for mdl in ['ols', 'rlm', 'rlm_pc5_x2']:
            cols = [zs(D.RTK), zs(l), zs(D.RTK) * zs(l), zs(D.purity), zs(D.prolif), zs(D.RTK) * zs(D.purity), zs(D.RTK) * zs(D.prolif)]
            if mdl == 'rlm_pc5_x2': cols += [zs(D.RTK)**2] + [zs(PCS.loc[D.index, f'PC{i}']) for i in range(1, 6)]
            M = sm.add_constant(np.column_stack(cols))
            r = sm.OLS(zs(D.ADAPT), M).fit() if mdl == 'ols' else sm.RLM(zs(D.ADAPT), M, M=sm.robust.norms.HuberT()).fit()
            rows.append(dict(gene=g, cohort=c, model=mdl, b=r.params[3], se=r.bse[3]))
R = pd.DataFrame(rows); o = []
for (g, mdl), d in R.groupby(['gene', 'model']):
    m = meta_re(d.b, d.se); loco = [meta_re(d[d.cohort != c].b, d[d.cohort != c].se)['p'] for c in d.cohort]
    o.append(dict(gene=g, model=mdl, est=m['est'], lo=m['est'] - 1.96 * m['se'], hi=m['est'] + 1.96 * m['se'], p=m['p'], I2=m['I2'], loco_max_p=max(loco)))
O = pd.DataFrame(o); O.to_csv(f'{OUT}/tables/T7_17q11_robustness.tsv', sep='\t', index=False); print(O.round(4).to_string())
