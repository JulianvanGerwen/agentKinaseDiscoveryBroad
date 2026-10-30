"""Controls: (a) split-half reliability, disattenuated and IV gains; (b) residual vs downstream-kinase protein (tc_padj and tc);
(c) spike-in power simulation for gain interactions on real data structure."""
from common import *
from scipy import stats
COV = ['purity', 'prolif']; rng = np.random.default_rng(3)
HALF = {'ADAPT': (S('SHC1', 'Y349', 'Y350', 'Y427'), S('PTPN11', 'Y542', 'Y580') + S('GAB1', 'Y627', 'Y659')),
        'ERK': (S('MAPK3', 'T202', 'Y204'), S('MAPK1', 'T185', 'Y187')),
        'RSK': (S('RPS6KA1', 'T359', 'S363', 'S380'), S('RPS6KA3', 'T365', 'S369', 'S386')),
        'AKTSUB': (S('AKT1S1', 'T246') + S('GSK3B', 'S9') + S('GSK3A', 'S21'), S('FOXO3', 'S253') + S('FOXO1', 'S256') + S('TSC2', 'S939', 'S981', 'T1462')),
        'MTORC1': (S('EIF4EBP1', 'T37', 'T46'), S('EIF4EBP1', 'S65', 'T70')),
        'RTK': (S('EGFR', 'Y1197', 'Y1110'), S('EGFR', 'Y1092', 'Y1172', 'Y1016') + S('ERBB2', 'Y1248', 'Y877') + S('MET', 'Y1234', 'Y1235')),
        'S6': (S('RPS6', 'S235', 'S236'), S('RPS6', 'S240', 'S244'))}
def resid(v, C):
    C1 = sm.add_constant(C); return v - C1 @ np.linalg.lstsq(C1, v, rcond=None)[0]
A = load_stages('tc_padj'); rows = []
for c in COHORTS:
    D = A[A.cohort == c]; X = cl.matrix(c, 'input_tc_padj')[D.index]
    halves = {}
    for st, (h1, h2) in HALF.items():
        a = X.reindex([s for s in h1 if s in X.index]).mean(0); b = X.reindex([s for s in h2 if s in X.index]).mean(0)
        halves[st] = (a.reindex(D.index).values.astype(float), b.reindex(D.index).values.astype(float))
    for e in EDGES:
        xs, ys = e
        if xs not in halves: continue
        a, b = halves[xs]; y = D[ys].values.astype(float); C = D[COV].values
        ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(y)
        if ok.sum() < 30: continue
        ra, rb, ry = resid(a[ok], C[ok]), resid(b[ok], C[ok]), resid(y[ok], C[ok])
        rab = np.corrcoef(ra, rb)[0, 1]; rel_x = 2 * rab / (1 + rab) if rab > 0 else np.nan
        rel_y = np.nan
        if ys in halves:
            ya, yb = halves[ys]; okk = ok & np.isfinite(ya) & np.isfinite(yb)
            if okk.sum() >= 30:
                q = np.corrcoef(resid(ya[okk], C[okk]), resid(yb[okk], C[okk]))[0, 1]; rel_y = 2 * q / (1 + q) if q > 0 else np.nan
        xc = resid(D[xs].values.astype(float)[ok], C[ok]); r_obs = np.corrcoef(xc, ry)[0, 1]
        # IV: y on half A instrumented by half B (and vice versa), standardized units of half A
        za, zb, zy = zs(ra), zs(rb), zs(ry)
        iv_ab = np.cov(zy, zb)[0, 1] / np.cov(za, zb)[0, 1]; iv_ba = np.cov(zy, za)[0, 1] / np.cov(zb, za)[0, 1]
        ols_a = np.corrcoef(za, zy)[0, 1]; ols_b = np.corrcoef(zb, zy)[0, 1]
        rows.append(dict(cohort=c, edge=ename(e), n=int(ok.sum()), half_r_x=rab, rel_x=rel_x, rel_y=rel_y, r_obs=r_obs,
                         r_disatt=r_obs / np.sqrt(rel_x * (rel_y if np.isfinite(rel_y) else 1)),
                         ols_halfA=ols_a, iv_halfA_by_B=iv_ab, ols_halfB=ols_b, iv_halfB_by_A=iv_ba))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/tables/T6_reliability_iv.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(R.groupby('edge')[['half_r_x', 'rel_x', 'rel_y', 'r_obs', 'r_disatt', 'ols_halfA', 'iv_halfA_by_B', 'ols_halfB', 'iv_halfB_by_A']].median().round(3).to_string())

# (b) residual vs downstream-kinase protein / total protein of output stage
RES = pd.read_parquet(f'{OUT}/tables/transmission_residuals.parquet')
PROT = {'ADAPT>ERK': ['P27361', 'P28482'], 'RTK>ERK': ['P27361', 'P28482'], 'ERK>RSK': ['Q15418', 'P51812'], 'ERK>ERKSUB': ['P27361', 'P28482'],
        'AKTSUB>MTORC1': ['Q13541', 'P42345'], 'MTORC1>S6': ['P23443', 'P62753'], 'ADAPT>AKTSUB': ['P31749', 'P31751'], 'RTK>ADAPT': ['P29353', 'Q06124']}
rows = []
for c in COHORTS:
    pr = cl.matrix(c, 'protein')
    for e, us in PROT.items():
        rr = RES[(RES.cohort == c) & (RES.edge == e)].set_index('sample_id')
        if len(rr) < 30: continue
        for u in us:
            if u not in pr.index: continue
            v = pr.loc[u].reindex(rr.index).values.astype(float); ok = np.isfinite(v)
            rows.append(dict(cohort=c, edge=e, protein=u, rho=stats.spearmanr(rr.resid.values[ok], v[ok])[0], n=int(ok.sum())))
P = pd.DataFrame(rows); P.to_csv(f'{OUT}/tables/T6_residual_vs_protein.tsv', sep='\t', index=False)
print(P.groupby(['edge', 'protein']).rho.agg(['median', 'min', 'max', 'size']).round(3).to_string())

# (c) spike-in power: y' = y + delta * x * L (standardized), L = real RNA of DUSP6 / PTEN; refit OLS interaction, meta
from common import meta_re
def int_fit(y, x, l, C):
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(l)
    y, x, l, C = zs(y[ok]), zs(x[ok]), zs(l[ok]), np.apply_along_axis(zs, 0, C[ok])
    M = sm.add_constant(np.column_stack([x, l, x * l, C, x[:, None] * C])); r = sm.OLS(y, M).fit(); return r.params[3], r.bse[3]
rows = []
for e, g in [(('ADAPT', 'ERK'), 'DUSP6'), (('ADAPT', 'AKTSUB'), 'PTEN'), (('AKTSUB', 'MTORC1'), 'TSC2')]:
    for delta in [0, 0.05, 0.1, 0.15, 0.2]:
        hits = 0; ests = []
        for rep in range(50):
            bs, ss = [], []
            for c in COHORTS:
                D = A[A.cohort == c]; rna = cl.matrix(c, 'rna')
                l = rna.loc[g].reindex(D.index).values.astype(float) if g in rna.index else None
                x = D[e[0]].values.astype(float); y = D[e[1]].values.astype(float)
                if l is None or np.isfinite(x * y).sum() < 30: continue
                ok = np.isfinite(x) & np.isfinite(y)
                yy = y.copy(); xs_ = (x - np.nanmean(x)) / np.nanstd(x); ls_ = (l - np.nanmean(l)) / np.nanstd(l)
                # perturb with random permutation of residual noise to make replicates differ
                yy[ok] = zs(y[ok]) + delta * xs_[ok] * ls_[ok] + rng.normal(0, 0.0, ok.sum())
                lperm = l if rep == 0 else l  # same L each rep; vary via bootstrap of tumours
                idx = rng.choice(np.where(ok)[0], ok.sum(), replace=True) if rep > 0 else np.where(ok)[0]
                b, s = int_fit(yy[idx], x[idx], l[idx], D[COV].values[idx]); bs.append(b); ss.append(s)
            m = meta_re(bs, ss); ests.append(m['est']); hits += m['p'] < 0.05
        rows.append(dict(edge=ename(e), gene=g, delta=delta, power=hits / 50, mean_est=np.mean(ests)))
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f'{OUT}/tables/T6_spikein_power.tsv', sep='\t', index=False)
