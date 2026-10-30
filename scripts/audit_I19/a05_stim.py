"""(f) Mirror test in non-inhibitor perturbation classes. Same metric as the inhibitor primary (experiment = study x cell x agent,
raw vs control, cross-study pairs, study-pair weighted). Prediction of 'direct portable / secondary context-specific':
activating perturbations (GF/RTK ligand, phosphatase inhibitors, DNA damage -> ATM/ATR) should show UP > DOWN replication.
Run on experiments and with time-matched conditions."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import *
OUT = '/workspace/results/audit_I19/'
pd.set_option('display.width', 250)
M = pd.read_csv(CA + 'stim_meta.tsv', sep='\t', dtype={'study_pmid': str})
drop = (M.control_raw == 'LPA') | (M.study_pmid == '29020638') | M.condition_raw.str.contains('p38i', na=False) | M.agent.isin(['cytarabine', 'methotrexate'])
M = M[~drop]
Xc = pd.read_parquet(CA + 'stim_cond.parquet')[M.conditionID]
M['exp'] = M.study_pmid + '|' + M.cell_line_norm.astype(str) + '|' + M.cls  # pool agents within class per study x cell (EGF/FGF etc.)
Xe = Xc.T.groupby(M.exp.values).mean().T
E = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), cls=('cls', 'first'), agent=('agent', lambda s: ','.join(sorted(set(map(str, s))))),
                         tmin=('time_min', 'min'), tmax=('time_min', 'max')).reindex(Xe.columns).reset_index()
keep = Xe.notna().sum(axis=1).values >= 3
A = Xe.values[keep].astype(np.float32)
# per-experiment global asymmetry (are stimulations up-skewed?)
E['frac_up1'] = np.nanmean(A > 1, 0); E['frac_dn1'] = np.nanmean(A < -1, 0); E['n_sites'] = (~np.isnan(A)).sum(0)
print(E.groupby('cls')[['frac_up1', 'frac_dn1']].median().round(3))
E.to_csv(OUT + 'T_stim_experiments.tsv', sep='\t', index=False)
rows = []
def run(label, R):
    if len(R) == 0: print(label, 'none'); return
    S = summarise(R)
    out = S[['dir', 'bin', 'pairs', 'study_pairs', 'rep', 'base', 'opp', 'opp_base', 'lor', 'meanB']].copy(); out.insert(0, 'test', label)
    out['dir_excess'] = (out.rep - out.base) - (out.opp - out.opp_base)
    for b in ['0.5-1', '1-2', '2-99']:
        try:
            o, lo, hi, ns, nsp = boot_diff(R, 'lor', b); out.loc[out.bin == b, 'dLOR'] = o; out.loc[out.bin == b, 'dLOR_CI'] = f'[{lo:.2f},{hi:.2f}]'
        except Exception: pass
    rows.append(out); print(out.round(3).to_string(), flush=True)
G = {c: (E.cls == c).values for c in E.cls.unique() if E[E.cls == c].study.nunique() >= 2}
print({c: (int(v.sum()), E[v].study.nunique()) for c, v in G.items()})
R = pair_rows(A, E, G, null_pool=None, min_both=100)
R.to_csv(OUT + 'cache/pairs_stim.tsv.gz', sep='\t', index=False)
act = ['GF_RTK', 'PPase_inh', 'DNA_damage', 'GPCR', 'TNF', 'PMA']
run('f_all_stim_classes', R)
run('f_activating(GF,PPase,DNA,GPCR,TNF,PMA)', R[R.group.isin(act)])
for c in G: run(f'f_{c}', R[R.group == c])
# early-only GF/RTK (<=60 min both)
early = ((E.tmax <= 60) | E.tmax.isna()).values
run('f_GF_RTK_early_only', pair_rows(A, E, {'GF_RTK': G['GF_RTK'] & early}, null_pool=None, min_both=100))
pd.concat(rows).to_csv(OUT + 'T_tests_f_stim.tsv', sep='\t', index=False)
