"""val_E step 2: re-estimate E abundance coefficient, IV decomposition, KRAS/EGFR offsets with ischaemia covariate."""
import importlib.util, sys, numpy as np, pandas as pd, warnings, statsmodels.api as sm; warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/audit_I24'); from common_a import *   # build, decomp, meta_re, zs, ACC, cl, COH
OUT = '/workspace/results/val_E/'
IS = pd.read_parquet('/workspace/cache/val_E/isch_scores.parquet')
def iscore(c, variant, col='isch_noPath'):
    d = IS[(IS.cohort == c) & (IS.variant == variant)]; return d[col]
def stage(c, v, isv='input_tc_padj'):
    D, _ = build(c, v); pr = cl.matrix(c, 'protein')
    for g in ['EGFR', 'ERBB2', 'MET']: D[g + '_prot'] = pr.loc[ACC[g]].reindex(D.index).values if ACC[g] in pr.index else np.nan
    D['RTK_prot'] = D[['EGFR_prot', 'ERBB2_prot', 'MET_prot']].mean(1)
    D['isch'] = iscore(c, isv).reindex(D.index); D['isch_raw'] = iscore(c, isv, 'isch').reindex(D.index)
    D['isch_tc'] = iscore(c, 'input_tc').reindex(D.index)
    pre = D.case_id.str.extract(r'^(\d+)')[0]; D['tss'] = pre.where(pre.notna() & (c in ['BRCA', 'COAD', 'OV']))
    D['cohort'] = c; return D
def ols_joint(D, covs, y='ERK', x=('RTK', 'RTK_prot'), tss=False):
    cols = [y, *x, *covs]; d = D[cols + (['tss'] if tss else [])].dropna()
    Z = d[cols].apply(zs); Xm = np.c_[np.ones(len(d)), Z[list(x) + list(covs)].values]
    if tss:
        T = pd.get_dummies(d.tss, drop_first=True).values.astype(float); Xm = np.c_[Xm, T]
    n, p = Xm.shape; b = np.linalg.lstsq(Xm, Z[y].values, rcond=None)[0]; r = Z[y].values - Xm @ b
    se = np.sqrt(np.diag(r @ r / (n - p) * np.linalg.inv(Xm.T @ Xm)))
    return dict(n=n, b_pY=b[1], se_pY=se[1], b_prot=b[2], se_prot=se[2])
ST = {(v, c): stage(c, v) for v in ['tc', 'tc_padj'] for c in COH}
rows = []
specs = {'base (purity, prolif)': ['purity', 'prolif'], '+isch (no pathway sites)': ['purity', 'prolif', 'isch'],
         '+isch raw (incl ERK1 T202)': ['purity', 'prolif', 'isch_raw'], '+isch from tc input': ['purity', 'prolif', 'isch_tc'],
         '+isch +isch^2': ['purity', 'prolif', 'isch', 'isch2']}
for c in COH:
    D = ST[('tc', c)].copy(); D['isch2'] = zs(D.isch) ** 2
    if D.RTK.notna().sum() < 30: continue
    for k, cv in specs.items():
        r = ols_joint(D, cv); r.update(cohort=c, spec=k); rows.append(r)
    if D.tss.notna().any():
        r = ols_joint(D, ['purity', 'prolif', 'isch'], tss=True); r.update(cohort=c, spec='+isch +site-ID proxy'); rows.append(r)
        r = ols_joint(D, ['purity', 'prolif'], tss=True); r.update(cohort=c, spec='base +site-ID proxy only'); rows.append(r)
J = pd.DataFrame(rows); J.to_csv(OUT + 'E_joint_by_cohort.tsv', sep='\t', index=False)
meta = []
for k, g in J.groupby('spec', sort=False):
    mp = meta_re(g.b_prot, g.se_prot); my = meta_re(g.b_pY, g.se_pY)
    meta.append(dict(spec=k, k=len(g), n=g.n.sum(), prot=mp['est'], prot_lo=mp['lo'], prot_hi=mp['hi'], prot_p=mp['p'], I2=mp['I2'],
                     n_neg=int((g.b_prot < 0).sum()), pY=my['est'], pY_lo=my['lo'], pY_hi=my['hi']))
Mj = pd.DataFrame(meta); Mj.to_csv(OUT + 'E_joint_meta.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(Mj.round(3).to_string())
print(J.pivot(index='cohort', columns='spec', values='b_prot').round(3).to_string())

# random-footprint negative control (over-adjustment): random 200-site score, same procedure, effect of adding on prot coefficient
rng = np.random.default_rng(1); ctl = []
Zc = {}
for c in COH:
    if ST[('tc', c)].RTK.notna().sum() < 30: continue
    x = cl.matrix(c, 'input_tc_padj'); x = x[[t for t in cl.tumour_ids(c) if t in x.columns]]; x = x[x.notna().mean(1) >= 0.5].astype(float)
    x = x[~x.index.str.split('_').str[0].isin(set(ACC.values()))]
    Zc[c] = x.sub(x.mean(1), axis=0).div(x.std(1), axis=0)
for it in range(30):
    rr = []
    for c, Z in Zc.items():
        D = ST[('tc', c)].copy(); idx = rng.choice(len(Z), 200, replace=False); sg = rng.choice([-1, 1], 200)
        D['rand'] = (Z.iloc[idx].mul(sg, axis=0)).mean(0).reindex(D.index)
        r = ols_joint(D, ['purity', 'prolif', 'rand']); r['cohort'] = c; rr.append(r)
    rr = pd.DataFrame(rr); m = meta_re(rr.b_prot, rr.se_prot); ctl.append(m['est'])
print('random-score control: meta RTK_prot coef mean %.3f sd %.3f range %.3f..%.3f' % (np.mean(ctl), np.std(ctl), min(ctl), max(ctl)))
pd.Series(ctl).to_csv(OUT + 'E_random_score_control.tsv', sep='\t')

# IV decomposition with isch covariate
rows = []
for v in ['tc', 'tc_padj']:
    for c in COH:
        D = ST[(v, c)]
        for name, cv in [('base', ('purity', 'prolif')), ('+isch', ('purity', 'prolif', 'isch'))]:
            r = decomp(D, cov=cv, nboot=500)
            if r: r.update(cohort=c, analysis=name, variant=v); rows.append(r)
R = pd.DataFrame(rows); R.to_csv(OUT + 'E_decomp_by_cohort.tsv', sep='\t', index=False)
M = meta_table(R, ['analysis', 'variant']); M.to_csv(OUT + 'E_decomp_meta.tsv', sep='\t', index=False)
print(M[['analysis', 'variant', 'k', 'medF', 'fs', 'rf', 'iv', 'iv_lo', 'iv_hi', 'other', 'other_lo', 'other_hi', 'diff', 'diff_p', 'total']].round(3).to_string())

# correlation of isch with receptor protein, RTK pY, KRAS
rows = []
for c in COH:
    D = ST[('tc', c)]; cm = cl.table(c, 'mutations'); mb = cl.matrix(c, 'mut_binary')
    kr = pd.Series(mb.loc['KRAS'].reindex(D.case_id).values, index=D.index) if 'KRAS' in mb.index else pd.Series(np.nan, index=D.index)
    kr[~D.case_id.isin(mb.columns)] = np.nan
    r = dict(cohort=c, n=len(D))
    for nme, o in [('EGFR_prot', D.EGFR_prot), ('ERBB2_prot', D.ERBB2_prot), ('RTK_prot', D.RTK_prot), ('RTK_pY', D.RTK), ('ADAPT_pY', D.ADAPT), ('ERK_pTpY', D.ERK), ('KRASmut', kr)]:
        dd = pd.concat([D.isch, o.astype(float)], axis=1).dropna()
        if len(dd) >= 20 and dd.iloc[:, 1].nunique() > 1:
            r[nme] = dd.corr(method='spearman').iloc[0, 1]; r['n_' + nme] = len(dd)
            if nme == 'KRASmut': r['nKRAS'] = int(dd.iloc[:, 1].sum())
    rows.append(r)
C = pd.DataFrame(rows); C.to_csv(OUT + 'isch_correlations.tsv', sep='\t', index=False)
print(C[['cohort', 'n', 'EGFR_prot', 'ERBB2_prot', 'RTK_prot', 'RTK_pY', 'ADAPT_pY', 'ERK_pTpY', 'KRASmut', 'nKRAS']].round(3).to_string())
print('median |rho|:', C[['EGFR_prot', 'ERBB2_prot', 'RTK_prot', 'RTK_pY', 'ADAPT_pY', 'ERK_pTpY', 'KRASmut']].abs().median().round(3).to_dict())
pd.to_pickle({k: v for k, v in ST.items() if k[0] == 'tc_padj'}, '/workspace/cache/val_E/st_padj.pkl')
