"""I29 shared regression helpers."""
import numpy as np, pandas as pd, statsmodels.api as sm, warnings
warnings.filterwarnings('ignore')
CA = '/workspace/cache/b_I29/'; OUT = '/workspace/results/b_I29/'
RNG = np.random.default_rng(29)

def z(s):
    s = pd.to_numeric(s, errors='coerce'); return (s - s.mean()) / s.std()

def load(coh='LUAD'):
    A = pd.read_parquet(CA + f'{coh}.parquet')
    A['ETVc'] = pd.concat([z(A['rna_' + g]) for g in ['ETV4', 'ETV5', 'SPRED1', 'DUSP6']], axis=1).mean(axis=1)
    A['pERK_broad'] = pd.concat([z(A.MAPK1_pTpY), z(A.MAPK3_pTpY)], axis=1).mean(axis=1)
    A['pERK_monoY'] = pd.concat([z(A.MAPK1_pY), z(A.MAPK3_pY)], axis=1).mean(axis=1)
    A['MAPKout_fp'] = -A['fp_KI_MEK']
    A['country'] = A.COUNTRY_OF_ORIGIN.where(A.COUNTRY_OF_ORIGIN.isin(['Vietnam', 'Usa', 'China']), 'Other').fillna('Other')
    A['log_excision'] = np.log(A.excision_min.clip(lower=1))
    A['log_excision_f'] = A.log_excision.fillna(A.log_excision.median())
    sm_ = A.SMOKING_STATUS.map({'Smoker': 1, 'Non-Smoker': 0})
    A['smoker'] = sm_.fillna(sm_.mean())
    st = A.STAGE.astype(str).str.extract(r'(IV|III|II|I)')[0].map({'I': 1, 'II': 2, 'III': 3, 'IV': 4})
    A['stage_num'] = st.fillna(st.median()) if st.notna().any() else 0
    A['male'] = (A.SEX.astype(str).str.lower().str.startswith('m')).astype(float)
    return A

def ols(A, y, covs=(), groups=('EGFR', 'RAS'), ref='WT', mask=None, fe=(), nperm=0, gcol='group'):
    d = A if mask is None else A[mask]
    d = d[d[gcol].isin(list(groups) + [ref])]
    cols = [y] + list(covs) + list(fe)
    d = d[[gcol] + cols].dropna()
    if len(d) < 10: return {}
    Y = z(d[y]).values  # z within analysed set
    Y = (d[y] - A[y].mean()) / A[y].std()  # z within cohort (all tumours)
    X = pd.DataFrame({g: (d[gcol] == g).astype(float) for g in groups}, index=d.index)
    for c in covs: X[c] = z(d[c]).fillna(0) if d[c].std() > 0 else 0
    for c in fe: X = X.join(pd.get_dummies(d[c].astype(str), prefix=c, drop_first=True).astype(float))
    X = X.loc[:, X.std() > 0] if len(X.columns) else X
    X = sm.add_constant(X, has_constant='add')
    r = sm.OLS(Y.values, X).fit(cov_type='HC3')
    out = {}
    for g in groups:
        if g not in X.columns: continue
        b, se = r.params[g], r.bse[g]
        o = dict(beta=b, lo=b - 1.96 * se, hi=b + 1.96 * se, p=r.pvalues[g], n_g=int((d[gcol] == g).sum()), n_ref=int((d[gcol] == ref).sum()), mde=2.8 * se)
        if nperm:
            # Freedman-Lane: permute residuals of reduced model
            Xr = X.drop(columns=[g]); rr = sm.OLS(Y.values, Xr).fit(); fit, res = rr.fittedvalues.values, rr.resid.values
            t0 = r.tvalues[g]; cnt = 0
            for _ in range(nperm):
                yp = fit + RNG.permutation(res)
                tp = sm.OLS(yp, X).fit(cov_type='HC3').tvalues[g]
                cnt += abs(tp) >= abs(t0)
            o['p_perm'] = (cnt + 1) / (nperm + 1)
        out[g] = o
    return out

def row(res, **kw):
    rows = []
    for g, o in res.items(): rows.append(dict(**kw, contrast=f'{g}_vs_ref', **o))
    return rows
