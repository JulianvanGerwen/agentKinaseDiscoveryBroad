"""Alternative-explanation tests (a-e) for the I19 up/down replication asymmetry. Primary variant: raw (vs control) values,
own-binder proteins masked, disjoint inhibitor classes with >=2 studies, cross-study experiment pairs, study-pair weighting."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import *
OUT = '/workspace/results/audit_I19/'
pd.set_option('display.width', 250)
X, E = load_exp('inh'); X = X[X.notna().sum(axis=1).values >= 6]
A = np.load(OUT + 'cache/A_raw.npy')
assert A.shape == X.shape
DM, HJ = '36926954', '31959955'
classes = [c for c in E.cls.unique() if E[E.cls == c].study.nunique() >= 2]
G = {c: (E.cls == c).values for c in classes}
allrows = []
def run(label, R, extra=''):
    S = summarise(R)
    out = S[['dir', 'bin', 'pairs', 'study_pairs', 'rep', 'base', 'snull', 'opp', 'opp_base', 'lor', 'meanB']].copy()
    out.insert(0, 'test', label)
    for b in ['1-2', '2-99']:
        try:
            o, lo, hi, ns, nsp = boot_diff(R, 'lor', b)
            out.loc[out.bin == b, 'dLOR_down_minus_up'] = o; out.loc[out.bin == b, 'dLOR_CI'] = f'[{lo:.2f},{hi:.2f}]'
        except Exception as e: pass
    allrows.append(out); print(out.round(3).to_string(), flush=True)
    return S
# T0 primary
R0 = pair_rows(A, E, G); run('T0_primary_raw_disjoint', R0)
# (a) detection / MNAR
nmeas_inh = (~np.isnan(A)).mean(1)
run('a1_sites_meas_>=50pct_of_310_exps', pair_rows(A, E, G, site_mask=nmeas_inh >= 0.5))
Rg = []
for c, idx in G.items():
    m = (~np.isnan(A[:, idx])).all(1)
    if m.sum() >= 200: Rg.append(pair_rows(A, E, {c: idx}, site_mask=m))
run('a2_sites_meas_in_ALL_exps_of_class', pd.concat(Rg))
Sn = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'n_conditions'], low_memory=False).set_index('site').n_conditions.reindex(X.index).values
run('a3_sites_n_conditions_top25pct', pair_rows(A, E, G, site_mask=Sn >= np.nanquantile(Sn, 0.75)))
# missingness in B conditional on discovery direction
mr = []
studies = E.study.values
for c, idx in G.items():
    ids = np.where(idx)[0]
    for i in ids:
        for j in ids:
            if studies[i] == studies[j]: continue
            a = A[:, i]; mj = np.isnan(A[:, j]); fin = ~np.isnan(a)
            mr.append(dict(si=studies[i], sj=studies[j], miss_all=mj[fin].mean(), miss_dn=mj[a < -1].mean() if (a < -1).sum() > 5 else np.nan,
                           miss_up=mj[a > 1].mean() if (a > 1).sum() > 5 else np.nan))
mr = pd.DataFrame(mr).groupby(['si', 'sj']).mean().mean()
print('missing in B | disc (study-pair mean):', mr.round(3).to_dict())
pd.DataFrame([mr]).to_csv(OUT + 'T_missingness.tsv', sep='\t', index=False)
# (d) composition
for lab, drop in [('d1_drop_decryptM', DM), ('d2_drop_Hijazi', HJ)]:
    keep = (E.study != drop).values
    run(lab, pair_rows(A, E, {c: idx & keep for c, idx in G.items()}))
c4 = [c for c in classes if E[E.cls == c].study.nunique() >= 4]
print('classes >=4 studies:', c4)
run('d3_classes_>=4_studies', R0[R0.group.isin(c4)])
for c in c4: run(f'd4_class_{c}', R0[R0.group == c])
# centred version for comparison
Ac = np.load(OUT + 'cache/A_centred.npy'); run('T0b_centred_disjoint', pair_rows(Ac, E, G))
pd.concat(allrows).to_csv(OUT + 'T_tests_a_d.tsv', sep='\t', index=False)
