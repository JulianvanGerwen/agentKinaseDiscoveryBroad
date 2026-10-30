"""H3 (+H4 exp-level): per-experiment axis shift scores vs reverse magnitude within study x cell units; cross-context concordance of shift."""
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, '/workspace/scripts/b_I25'); from common_i25 import *
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import load_cond
OUT = '/workspace/results/b_I25/tables/'
rng = np.random.default_rng(3)
def shifts(A, AX):
    """A sites x exps -> DataFrame exps x features."""
    V = A.values; out = {}
    def mmean(mask):
        M = V[mask]; n = np.sum(~np.isnan(M), 0); return np.where(n >= 10, np.nanmean(M, 0), np.nan)
    out['shift_M'] = mmean(AX.CC_Mup.values) - mmean(AX.CC_Mdown.values)
    for ax, sc in [('MKI', 'MKI_score'), ('INTER', 'INTER_score'), ('MET24', 'MET24_score'), ('MET30', 'MET30_score'), ('CULT3D', 'CULT3D_score'), ('SERUM', 'SERUM_score'), ('OIS', 'OIS_score')]:
        s = AX[sc].values; ok = (np.abs(np.nan_to_num(s)) >= 1.5)
        if ax == 'MKI': out['shift_MKI'] = mmean(ok & (s > 0)) - mmean(AX.U_MKI.values & ~ok)
        else: out[f'shift_{ax}'] = mmean(ok & (s > 0)) - mmean(ok & (s < 0))
    for st in ['S_HIPPO', 'S_AMPKmotif']:
        out['shift_' + st[2:]] = mmean(AX[st].values) - mmean(~AX[st].values)
    nc = ~AX.CC_union.values; Vn = V[nc]; meas = np.sum(~np.isnan(Vn), 0)
    out['R_rev'] = np.nansum(Vn > 1, 0) / meas; out['P_prim'] = np.nansum(Vn < -1, 0) / meas
    out['frac_down'] = np.nansum(V < -0.3, 0) / np.sum(~np.isnan(V), 0); out['median'] = np.nanmedian(V, 0)
    out['sd'] = np.nanstd(V, 0)
    return pd.DataFrame(out, index=A.columns)
A, E = load_test('inh'); AX = axes(A.index)
F = shifts(A, AX); F = pd.concat([E.set_index('exp'), F], axis=1); F.index.name = 'exp'
F.to_csv(OUT + 'T3_H3_exp_features.tsv', sep='\t')
FEATS = [c for c in F.columns if c.startswith('shift_')] + ['frac_down', 'median']
def partial_rho(d, x, y='R_rev', cov=('P_prim', 'frac_down')):
    d = d[[x, y, *cov]].dropna()
    if len(d) < 6: return np.nan, np.nan, len(d)
    r = d.rank(); Z = np.c_[np.ones(len(r)), r[list(cov)].values]
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    cov_eff = [c for c in cov if c != x]
    Z = np.c_[np.ones(len(r)), r[cov_eff].values] if cov_eff else np.ones((len(r), 1))
    return stats.spearmanr(d[x], d[y])[0], np.corrcoef(res(r[x].values), res(r[y].values))[0, 1], len(d)
F['unit'] = F.study + '|' + F.cell
units = [u for u, d in F.groupby('unit') if len(d) >= 6]
rows = []
for f in FEATS:
    per = []
    for u in units:
        rho, prho, n = partial_rho(F[F.unit == u], f, cov=('P_prim', 'frac_down') if f != 'frac_down' else ('P_prim',))
        per.append((u, rho, prho, n))
    P = pd.DataFrame(per, columns=['unit', 'rho', 'prho', 'n']).dropna()
    # absolute shift version
    pa = []
    for u in units:
        d = F[F.unit == u].copy(); d['abs'] = d[f].abs(); pa.append(partial_rho(d, 'abs')[1])
    bs = [P.sample(len(P), replace=True, random_state=k)[['rho', 'prho']].mean().values for k in range(1000)]
    bs = np.array(bs)
    rows.append(dict(feature=f, n_units=len(P), n_exp=int(P.n.sum()), mean_rho=P.rho.mean(), mean_prho=P.prho.mean(),
                     prho_lo=np.percentile(bs[:, 1], 2.5), prho_hi=np.percentile(bs[:, 1], 97.5), rho_lo=np.percentile(bs[:, 0], 2.5), rho_hi=np.percentile(bs[:, 0], 97.5),
                     mean_prho_abs=np.nanmean(pa), units_pos=int((P.prho > 0).sum())))
    for r in per: rows_u = None
T3 = pd.DataFrame(rows); T3.to_csv(OUT + 'T3_H3_within_unit_corr.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(T3.round(3).to_string())
# per-unit table for shift_M
pu = []
for u in units:
    d = F[F.unit == u]; r = partial_rho(d, 'shift_M'); pu.append(dict(unit=u, n=r[2], rho=r[0], prho=r[1], shift_M_mean=d.shift_M.mean(), shift_M_sd=d.shift_M.std()))
print(pd.DataFrame(pu).round(3).to_string())
pd.DataFrame(pu).to_csv(OUT + 'T3b_H3_shiftM_per_unit.tsv', sep='\t', index=False)
# ---- concordance of shifts across contexts
# (a) within-experiment replicate concordance: split conditions of the same experiment (different dose/time)
Xc, Ec = load_cond('inh'); Ec = Ec[Ec.exp.isin(F.index)]; Xc = Xc[Ec.conditionID]; Xc = Xc.reindex(A.index)
Fc = shifts(Xc, AX); Fc['exp'] = Ec.exp.values
conc = []
def pairs_corr(feature):
    out = {}
    g = Fc.groupby('exp')[feature]
    a = g.nth(0); b = g.nth(1); ex = Fc.groupby('exp').size(); ex = ex[ex >= 2].index
    a = Fc.groupby('exp')[feature].apply(lambda s: s.iloc[0]).reindex(ex); b = Fc.groupby('exp')[feature].apply(lambda s: s.iloc[-1]).reindex(ex)
    m = a.notna() & b.notna(); out['within_exp_r'] = np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 5 else np.nan; out['within_exp_n'] = int(m.sum())
    # (b) same drug, different cell, same study (Hijazi)
    H = F[F.study == '31959955'].pivot_table(index='agent', columns='cell', values=feature)
    rs = []
    for i, c1 in enumerate(H.columns):
        for c2 in H.columns[i + 1:]:
            m = H[c1].notna() & H[c2].notna(); rs.append(np.corrcoef(H.loc[m, c1], H.loc[m, c2])[0, 1])
    out['hijazi_samedrug_diffcell_r'] = np.mean(rs)
    # (c) same class, different study: class means per unit, correlate across unit pairs over shared classes (>=4)
    Q = F.pivot_table(index='cls', columns='unit', values=feature); rs = []
    for i, u1 in enumerate(Q.columns):
        for u2 in Q.columns[i + 1:]:
            if u1.split('|')[0] == u2.split('|')[0]: continue
            m = Q[u1].notna() & Q[u2].notna()
            if m.sum() >= 4: rs.append(np.corrcoef(Q.loc[m, u1], Q.loc[m, u2])[0, 1])
    out['xstudy_sameclass_r'] = np.nanmean(rs) if rs else np.nan; out['xstudy_n_unitpairs'] = len(rs)
    # (d) variance decomposition: unit vs class (two-way OLS, type-I eta2 both orders averaged)
    import statsmodels.formula.api as smf
    d = F[[feature, 'unit', 'cls']].dropna().rename(columns={feature: 'y'})
    tot = ((d.y - d.y.mean()) ** 2).sum()
    r_u = smf.ols('y ~ C(unit)', d).fit().ess / tot; r_c = smf.ols('y ~ C(cls)', d).fit().ess / tot; r_b = smf.ols('y ~ C(unit)+C(cls)', d).fit().ess / tot
    out.update(R2_unit=r_u, R2_class=r_c, R2_both=r_b, uniq_unit=r_b - r_c, uniq_class=r_b - r_u)
    return out
for f in ['shift_M', 'shift_MKI', 'shift_INTER', 'shift_MET24', 'shift_CULT3D', 'shift_SERUM', 'shift_OIS', 'shift_HIPPO', 'shift_AMPKmotif', 'frac_down', 'P_prim', 'R_rev']:
    r = pairs_corr(f); r['feature'] = f; conc.append(r)
T3c = pd.DataFrame(conc); T3c.to_csv(OUT + 'T3c_shift_concordance.tsv', sep='\t', index=False)
print(T3c.round(3).to_string())
