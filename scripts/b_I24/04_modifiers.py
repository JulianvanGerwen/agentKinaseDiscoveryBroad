"""Pre-declared gain modifiers: y ~ x + L + x:L + purity + prolif + x:purity + x:prolif (robust, standardized), per cohort; RE meta.
Sensitivity model adds x^2 and 5 global phospho PCs."""
from common import *
from statsmodels.stats.multitest import multipletests
UNI = dict(DUSP4='Q13115', DUSP6='Q16828', DUSP5='Q16690', SPRY2='O43597', SPRY4='Q9C004', SPRED1='Q7Z699', SPRED2='Q7Z698',
           NF1='P21359', RASA1='P20936', PTPN11='Q06124', PTEN='P60484', INPP4B='O15327', PHLPP1='O60346', PHLPP2='Q6ZVD8',
           TSC1='Q92574', TSC2='P49815', STK11='Q15831', PPP2R2A='P63151', KRAS='P01116', BRAF='P15056', EGFR='P00533', PIK3CA='P42336')
MAPKREG = ['DUSP4', 'DUSP6', 'DUSP5', 'SPRY2', 'SPRY4', 'SPRED1', 'SPRED2', 'NF1', 'RASA1', 'PTPN11', 'KRAS', 'BRAF', 'EGFR', 'PPP2R2A']
PRIMARY = {'ADAPT>ERK': MAPKREG, 'RTK>ERK': MAPKREG, 'ERK>RSK': ['DUSP4', 'DUSP5', 'DUSP6', 'PPP2R2A'], 'ERK>ERKSUB': ['DUSP4', 'DUSP5', 'DUSP6', 'PPP2R2A'],
           'RTK>ADAPT': ['EGFR', 'PTPN11'], 'ADAPT>AKTSUB': ['PTEN', 'INPP4B', 'PHLPP1', 'PHLPP2', 'PIK3CA', 'PPP2R2A'],
           'AKTSUB>MTORC1': ['TSC1', 'TSC2', 'STK11']}
# predicted sign of interaction for continuous abundance (negative regulators: more -> lower gain); mutation = LoF for suppressors
NEG = {'DUSP4', 'DUSP6', 'DUSP5', 'SPRY2', 'SPRY4', 'SPRED1', 'SPRED2', 'NF1', 'RASA1', 'PTEN', 'INPP4B', 'PHLPP1', 'PHLPP2', 'TSC1', 'TSC2', 'STK11', 'PPP2R2A'}
COV = ['purity', 'prolif']
A = load_stages('tc_padj'); PCS = pd.read_parquet(f'{OUT}/tables/phospho_pcs.parquet')

def lesions(c, D):
    cases = D.case_id.values; L = {}
    cn = cl.matrix(c, 'cnv_log2'); rna = cl.matrix(c, 'rna'); pr = cl.matrix(c, 'protein'); mb = cl.matrix(c, 'mut_binary'); gi = cl.matrix(c, 'cnv_gistic')
    for g, u in UNI.items():
        if g in cn.index: L[(g, 'cnv')] = cn.loc[g].reindex(cases).values.astype(float)
        if g in rna.index: L[(g, 'rna')] = rna.loc[g].reindex(D.index).values.astype(float)
        if u in pr.index: L[(g, 'protein')] = pr.loc[u].reindex(D.index).values.astype(float)
        mcases = [x for x in cases if x in mb.columns]
        v = (mb.loc[g].reindex(cases).fillna(0).values if g in mb.index else np.zeros(len(cases))).astype(float)
        v[~np.isin(cases, mcases)] = np.nan
        if np.nansum(v) >= 5: L[(g, 'mut')] = v
    if 'EGFR' in gi.index:
        v = (gi.loc['EGFR'].reindex(cases).values == 2).astype(float)
        if v.sum() >= 5: L[('EGFR', 'amp')] = v
    return L

def fit(y, x, l, C, binary, extra=None):
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(l) & np.all(np.isfinite(C), 1)
    if ok.sum() < 30: return None
    if binary and (l[ok].sum() < 5 or (1 - l[ok]).sum() < 5): return None
    y_, x_, l_, C_ = zs(y[ok]), zs(x[ok]), l[ok] if binary else zs(l[ok]), np.apply_along_axis(zs, 0, C[ok])
    cols = [x_, l_, x_ * l_] + [C_[:, j] for j in range(C_.shape[1])] + [x_ * C_[:, j] for j in range(min(2, C_.shape[1]))]
    if extra is not None:
        E = np.apply_along_axis(zs, 0, extra[ok]); cols += [x_**2] + [E[:, j] for j in range(E.shape[1])]
    Xm = sm.add_constant(np.column_stack(cols))
    try: r = sm.RLM(y_, Xm, M=sm.robust.norms.HuberT()).fit()
    except Exception: return None
    return dict(n=int(ok.sum()), n_mut=int(l[ok].sum()) if binary else np.nan, b_int=r.params[3], se_int=r.bse[3], b_lev=r.params[2], se_lev=r.bse[2],
                b_x=r.params[1])

rows = []
for c in COHORTS:
    D = A[A.cohort == c]; L = lesions(c, D); C = D[COV].values; P5 = PCS.loc[D.index].values[:, :5]
    for e in EDGES:
        xs, ys = e; y, x = D[ys].values.astype(float), D[xs].values.astype(float)
        if np.isfinite(y).sum() < 30 or np.isfinite(x).sum() < 30: continue
        for (g, t), l in L.items():
            b = t in ('mut', 'amp')
            for model, extra in [('main', None), ('sens_pc_x2', P5)]:
                r = fit(y, x, l, C, b, extra)
                if r: rows.append(dict(cohort=c, edge=ename(e), gene=g, type=t, model=model, **r))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/tables/T4_modifiers_by_cohort.tsv', sep='\t', index=False)
out = []
for (e, g, t, mdl), d in R.groupby(['edge', 'gene', 'type', 'model']):
    mi = meta_re(d.b_int, d.se_int); ml = meta_re(d.b_lev, d.se_lev)
    zc = d.b_int / d.se_int
    out.append(dict(edge=e, gene=g, type=t, model=mdl, k=mi['k'], int_meta=mi['est'], int_lo=mi['est'] - 1.96 * mi['se'], int_hi=mi['est'] + 1.96 * mi['se'],
                    int_p=mi['p'], int_I2=mi['I2'], n_nom_same_sign=int(((np.sign(zc) == np.sign(mi['est'])) & (abs(zc) > 1.96)).sum()),
                    lev_meta=ml['est'], lev_lo=ml['est'] - 1.96 * ml['se'], lev_hi=ml['est'] + 1.96 * ml['se'], lev_p=ml['p'],
                    primary=g in PRIMARY.get(e, []),
                    predicted_sign=(-1 if g in NEG else 1) * (1 if t in ('cnv', 'rna', 'protein') else -1) if g in NEG else (1 if t in ('mut', 'amp') else np.nan)))
M = pd.DataFrame(out)
for mdl in ['main', 'sens_pc_x2']:
    m = (M.model == mdl)
    M.loc[m, 'int_q_all'] = multipletests(M.loc[m, 'int_p'].fillna(1), method='fdr_bh')[1]
    M.loc[m, 'lev_q_all'] = multipletests(M.loc[m, 'lev_p'].fillna(1), method='fdr_bh')[1]
    mp = m & M.primary
    M.loc[mp, 'int_q_primary'] = multipletests(M.loc[mp, 'int_p'].fillna(1), method='fdr_bh')[1]
    M.loc[mp, 'lev_q_primary'] = multipletests(M.loc[mp, 'lev_p'].fillna(1), method='fdr_bh')[1]
M.to_csv(f'{OUT}/tables/T4_modifiers_meta.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
cols = ['edge', 'gene', 'type', 'k', 'int_meta', 'int_lo', 'int_hi', 'int_p', 'int_q_primary', 'n_nom_same_sign', 'int_I2', 'lev_meta', 'lev_p', 'predicted_sign']
print('PRIMARY main, sorted by int p'); print(M[(M.model == 'main') & M.primary].sort_values('int_p')[cols].head(30).round(4).to_string())
print('ALL main top int'); print(M[(M.model == 'main')].sort_values('int_p')[cols + ['int_q_all']].head(25).round(4).to_string())
print('n primary tests', ((M.model == 'main') & M.primary).sum(), 'n all', (M.model == 'main').sum())
