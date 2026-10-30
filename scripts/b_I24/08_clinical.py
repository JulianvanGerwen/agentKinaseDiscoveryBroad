"""Transmission residuals vs subtype (Kruskal), therapy-relevant states, and OS (Cox PHReg, per cohort, RE meta)."""
from common import *
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.duration.hazard_regression import PHReg
import os
RES = pd.read_parquet(f'{OUT}/tables/transmission_residuals.parquet')
EDG = ['RTK>ADAPT', 'ADAPT>ERK', 'ERK>RSK', 'ERK>ERKSUB', 'ADAPT>AKTSUB', 'AKTSUB>MTORC1', 'MTORC1>S6']
W = RES[RES.edge.isin(EDG)].pivot(index='sample_id', columns='edge', values='resid')
W['cohort'] = RES.drop_duplicates('sample_id').set_index('sample_id').cohort.reindex(W.index); W['case_id'] = [s[:-2] for s in W.index]
SUBCOL = {'BRCA': 'PAM50', 'LUAD': 'NMF_CONSENSUS', 'LSCC': 'NMF_SUBTYPE', 'COAD': 'MSI_STATUS', 'UCEC': 'HISTOLOGIC_TYPE'}
rows = []
for c, col in SUBCOL.items():
    f = f'/workspace/results/infra_cptac/cohorts/{c}/subtypes_cbioportal.parquet'
    if not os.path.exists(f): continue
    s = pd.read_parquet(f).drop_duplicates('case_id').set_index('case_id')
    if col not in s.columns: continue
    d = W[W.cohort == c].copy(); d['sub'] = s[col].reindex(d.case_id).values
    for e in EDG:
        q = d[[e, 'sub']].dropna(); vc = q['sub'].value_counts(); keep = vc[vc >= 5].index; q = q[q['sub'].isin(keep)]
        if q['sub'].nunique() < 2: continue
        H, p = stats.kruskal(*[g[e].values for _, g in q.groupby('sub')])
        means = q.groupby('sub')[e].mean().round(2).to_dict()
        rows.append(dict(cohort=c, subtype=col, edge=e, n=len(q), H=H, p=p, group_means=str(means)))
S_ = pd.DataFrame(rows); S_['q'] = multipletests(S_.p, method='fdr_bh')[1]
S_.to_csv(f'{OUT}/tables/T8_subtype.tsv', sep='\t', index=False); pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 120)
print(S_.sort_values('p').head(12).to_string())
# survival
clin = cl.clinical().drop_duplicates('case_id').set_index('case_id'); rows = []
for c in COHORTS:
    d = W[W.cohort == c].copy(); d['t'] = clin.OS_days.reindex(d.case_id).values; d['ev'] = clin.OS_event.reindex(d.case_id).values
    d['age'] = clin.Age.reindex(d.case_id).values
    for e in EDG:
        q = d[[e, 't', 'ev']].dropna(); q = q[q.t > 0]
        if len(q) < 30 or q.ev.sum() < 8: continue
        m = PHReg(q.t.values, zs(q[e].values)[:, None], status=q.ev.values.astype(int)).fit()
        rows.append(dict(cohort=c, edge=e, n=len(q), events=int(q.ev.sum()), logHR=m.params[0], se=m.bse[0]))
SV = pd.DataFrame(rows); SV.to_csv(f'{OUT}/tables/T8_survival_by_cohort.tsv', sep='\t', index=False)
o = [dict(edge=e, **meta_re(g.logHR, g.se)) for e, g in SV.groupby('edge')]
O = pd.DataFrame(o); O['HR'] = np.exp(O.est); O['q'] = multipletests(O.p, method='fdr_bh')[1]
O.to_csv(f'{OUT}/tables/T8_survival_meta.tsv', sep='\t', index=False); print(O.round(3).to_string())
