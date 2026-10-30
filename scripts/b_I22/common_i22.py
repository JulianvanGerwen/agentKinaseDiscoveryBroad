import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import load_exp, mask_binders, CA as ACA
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'residue', 'gene'], low_memory=False).set_index('site')

def get(kind, st_only=True):
    if kind == 'inh':
        X, E = load_exp('inh'); A = pd.DataFrame(mask_binders(X, E), index=X.index, columns=X.columns)
        E = E.rename(columns={E.columns[0]: 'exp'}); E = E[['exp', 'study', 'cell', 'cls', 'agent', 'time_min', 'dose_nM']].copy()
    else:
        M = pd.read_csv(ACA + 'stim_meta.tsv', sep='\t', dtype={'study_pmid': str})
        drop = (M.control_raw == 'LPA') | (M.study_pmid == '29020638') | M.condition_raw.str.contains('p38i', na=False) | M.agent.isin(['cytarabine', 'methotrexate'])
        M = M[~drop]; Xc = pd.read_parquet(ACA + 'stim_cond.parquet')[M.conditionID]
        M['exp'] = M.study_pmid + '|' + M.cell_line_norm.astype(str) + '|' + M.cls
        A = Xc.T.groupby(M.exp.values).mean().T
        E = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), cls=('cls', 'first'), agent=('agent', 'first'),
                                 time_min=('time_min', 'median')).reindex(A.columns).reset_index()
        E = E.rename(columns={E.columns[0]: 'exp'}); E['dose_nM'] = np.nan
    if st_only: A = A[S.residue.reindex(A.index).isin(['S', 'T']).values]
    E['kind'] = kind; E['cell'] = E.cell.astype(str)
    return A, E.reset_index(drop=True)
