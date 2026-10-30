"""Split unadjusted RTK pY into receptor-abundance component (fitted from RTK protein) and stoichiometry component
(residual); IV-decompose adaptor->ERK through each. Also AKT-substrate outcome."""
from common_a import *
pd.set_option('display.width', 250)
A = pd.read_parquet(f'{OUT}/tables/stages_audit.parquet')
AKT = ['Q96B36_T246', 'P49841_S9', 'P49840_S21', 'O43524_S253', 'Q12778_S256', 'P49815_S939', 'P49815_S981', 'P49815_T1462']
rows = []
for v in ['tc', 'tc_padj']:
    for c in COH:
        D = A[(A.variant == v) & (A.cohort == c)].copy()
        X = cl.matrix(c, f'input_{v}')[D.index]; D['AKTSUB'] = X.reindex([s for s in AKT if s in X.index]).mean(0)
        if D.RTK.notna().mean() < 0.4: continue
        r = decomp(D, y='AKTSUB')
        if r: r.update(cohort=c, analysis='AKTSUB_outcome', variant=v); rows.append(r)
        if v == 'tc':
            d = D[D.RTK.notna() & D.RTK_prot.notna()].copy()
            b = np.polyfit(d.RTK_prot, d.RTK, 1); d['RTK_abund'] = np.polyval(b, d.RTK_prot); d['RTK_stoich'] = d.RTK - d.RTK_abund
            for x in ['RTK_abund', 'RTK_stoich']:
                r = decomp(d, x=x)
                if r: r.update(cohort=c, analysis=f'instr_{x}', variant=v); rows.append(r)
            # reduced-form joint: ERK ~ RTKpY + RTKprot
            dd = d[['ERK', 'RTK', 'RTK_prot', 'ADAPT', 'purity', 'prolif']].dropna().apply(zs)
            for y in ['ERK', 'ADAPT']:
                Xm = np.c_[np.ones(len(dd)), dd[['RTK', 'RTK_prot', 'purity', 'prolif']].values]
                bb = np.linalg.lstsq(Xm, dd[y].values, rcond=None)[0]; res = dd[y].values - Xm @ bb
                se = np.sqrt(np.diag(res @ res / (len(dd) - 5) * np.linalg.inv(Xm.T @ Xm)))
                rows.append(dict(cohort=c, analysis=f'joint_{y}~RTKpY+RTKprot', variant=v, n=len(dd), fs=bb[1], fs_se=se[1], rf=bb[2], rf_se=se[2]))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/tables/A3_abundance_split_by_cohort.tsv', sep='\t', index=False)
M = meta_table(R[R.analysis.str.startswith(('AKTSUB', 'instr'))], ['analysis', 'variant'])
M.to_csv(f'{OUT}/tables/A3_abundance_split_meta.tsv', sep='\t', index=False)
print(M[['analysis', 'variant', 'k', 'medF', 'minF', 'fs', 'rf', 'iv', 'iv_lo', 'iv_hi', 'other', 'other_lo', 'other_hi', 'diff', 'diff_p']].round(3).to_string())
print(R[R.analysis.str.startswith('instr')][['analysis', 'cohort', 'n', 'F', 'fs', 'rf', 'iv', 'iv_se', 'other']].round(3).to_string())
J = R[R.analysis.str.startswith('joint')]
for a, g in J.groupby('analysis'):
    m1, m2 = meta_re(g.fs, g.fs_se), meta_re(g.rf, g.rf_se)
    print(f'{a}: coef RTK pY {m1["est"]:.3f} [{m1["lo"]:.3f},{m1["hi"]:.3f}] ; coef RTK protein {m2["est"]:.3f} [{m2["lo"]:.3f},{m2["hi"]:.3f}] p={m2["p"]:.3g}')
    print(g[['cohort', 'fs', 'rf']].round(3).to_string())
