"""Task 2: portability vs timing, with the audit_I19 cross-study replication metric.
Discovery = one time series (study S). Replication = audit_I19 experiments (study x cell x agent; stimulations pooled per class
as in a05) of the matching perturbation class from studies != S (no condition is used both to define and to test).
Site categories in the discovery series: direction (primary / reverse, sign of peak), timing (onset = first time point reaching
half of the peak; 'early' = first sampled time point, 'late' otherwise; early split into sustained (still > half-peak... at last tp)
/ transient), |peak| bin. Replication in exp j: dir*x_j > 0.3; LOR vs non-discovery co-measured sites (per-direction base rate),
exactly the audit pair metric. Summaries: mean over pairs within study pair, then over study pairs; CIs by study bootstrap."""
import sys, pickle, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import load_exp, mask_binders, CA as ACA
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
SM = pd.read_csv(CA + 'series_meta.tsv', sep='\t', dtype={'study': str})
SM = SM[SM.prim.isin(['up', 'down']) & SM.rcls.notna()]
L = pd.read_parquet(CA + 'onset_long.parquet')
VV = pickle.load(open(CA + 'series_V.pkl', 'rb'))
# replication experiments
Xi, Ei = load_exp('inh'); Ai = pd.DataFrame(mask_binders(Xi, Ei), index=Xi.index); Ei['kind'] = 'inh'
M = pd.read_csv(ACA + 'stim_meta.tsv', sep='\t', dtype={'study_pmid': str})
drop = (M.control_raw == 'LPA') | (M.study_pmid == '29020638') | M.condition_raw.str.contains('p38i', na=False) | M.agent.isin(['cytarabine', 'methotrexate'])
M = M[~drop]; Xc = pd.read_parquet(ACA + 'stim_cond.parquet')[M.conditionID]
M['exp'] = M.study_pmid + '|' + M.cell_line_norm.astype(str) + '|' + M.cls
As = Xc.T.groupby(M.exp.values).mean().T
Es = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), cls=('cls', 'first')).reindex(As.columns).reset_index()
Es['kind'] = 'stim'
EXP = {'inh': (Ai, Ei), 'stim': (As, Es)}
BIN = lambda a: np.select([a <= 1, a <= 2], ['0.5-1', '1-2'], '2-99')
rows = []
for s in SM.itertuples():
    kind = 'inh' if s.prim == 'down' else 'stim'
    A, E = EXP[kind]
    js = np.where((E.cls == s.rcls).values & (E.study != s.study).values)[0]
    if len(js) == 0: continue
    sites, V, times = VV[s.series]
    pk = V[np.arange(len(V)), np.abs(V).argmax(1)]; dv = V * np.sign(pk)[:, None]
    on_half = (dv >= 0.5 * np.abs(pk)[:, None]).argmax(1)
    on1 = np.where((dv > 1).any(1), (dv > 1).argmax(1), -1)
    sus = dv[:, -1] >= 0.5 * np.abs(pk)
    sg = 1 if s.prim == 'up' else -1
    D = pd.DataFrame(dict(site=sites, peak=pk, dirv=np.sign(pk), prim=np.sign(pk) == sg, on_half=on_half, on1=on1, sus=sus))
    D = D[D.peak.abs() > 0.5]
    D['bin'] = BIN(D.peak.abs().values)
    D['timing'] = np.where(D.on_half > 0, 'late', np.where(D.sus, 'early_sus', 'early_trans'))
    D['timing1'] = np.where(D.peak.abs() <= 1, 'na', np.where(D.on1 > 0, 'late', 'early'))
    sidx = pd.Index(sites)
    for j in js:
        xj = A.iloc[:, j].reindex(sites).values
        fin = ~np.isnan(xj)
        if fin.sum() < 200: continue
        base_rep = {d: (d * xj[fin] > 0.3) for d in (1, -1)}
        Df = D[fin[sidx.get_indexer(D.site)]]
        xjD = xj[sidx.get_indexer(Df.site)]
        for tcol in ['timing', 'timing1']:
            for (pr, tm, b), g in Df.groupby(['prim', tcol, 'bin']):
                if tm == 'na': continue
                for d, gg in g.groupby('dirv'):
                    nd = len(gg)
                    if nd < 5: continue
                    y = d * xjD[Df.index.get_indexer(gg.index)] > 0.3
                    allrep = base_rep[d]; tot = fin.sum()
                    a_ = y.sum() + .5; b_ = nd - y.sum() + .5; c_ = allrep.sum() - y.sum() + .5; d_ = (tot - nd) - (allrep.sum() - y.sum()) + .5
                    rows.append(dict(series=s.series, si=s.study, sj=E.study.iloc[j], cell_i=s.cell, cell_j=E.cell.iloc[j], same_cell=str(s.cell) == str(E.cell.iloc[j]),
                                     rcls=s.rcls, kind=kind, tdef=tcol, prim='primary' if pr else 'reverse', timing=tm, bin=b, dirv=int(d), n_disc=nd,
                                     rep=y.mean(), base=allrep.mean(), lor=np.log(a_ * d_ / (b_ * c_))))
R = pd.DataFrame(rows)
R.to_csv(CA + 'portability_pairs.tsv.gz', sep='\t', index=False)
print(R.shape, R.groupby('kind').si.nunique(), R.groupby('kind').sj.nunique(), R.same_cell.mean())
