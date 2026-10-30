"""PLACEBO for 07 (Task 3a): class kinase set replaced by random same-size kinase sets (5 draws per class).
Original doc: among primary-direction sites, are target-motif-matched (putatively direct) sites more portable?
All audit_I19 experiments (inhibitors: binder-masked; stimulations pooled per class), cross-study ordered pairs in the same class.
Discovery in exp i (dir*x_i in |log2FC| bin) split into motif-matched (max Cantley pct over the class's direct kinases >= 90,
residue-class appropriate) vs not; replication in j (dir*x_j > 0.3), LOR vs all non-discovery co-measured sites (same baseline).
Also PSP known substrates of the class kinases. Reverse direction = negative control for the motif split (and the question
whether reverse sites matching the target motif behave like primary ones).
Motif is sequence-only -> independent of both i and j."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import load_exp, mask_binders, CA as ACA
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
KMAP = {'MEK': 'MAPK1 MAPK3', 'ERK': 'MAPK1 MAPK3', 'RAF': 'MAPK1 MAPK3', 'EGFR': 'EGFR MAPK1 MAPK3', 'PI3K': 'AKT1 AKT2 AKT3',
        'PI3K_mTOR': 'AKT1 MTOR RPS6KB1', 'mTORcat': 'MTOR RPS6KB1 AKT1', 'rapalog': 'MTOR RPS6KB1', 'AKT': 'AKT1 AKT2 AKT3', 'PDPK1': 'PDPK1 AKT1 RPS6KB1',
        'S6K': 'RPS6KB1 RPS6KB2', 'SRC_ABL': 'SRC ABL1', 'FGFR': 'FGFR1 FGFR2 MAPK1 MAPK3', 'MET_ALK': 'MET ALK MAPK1 MAPK3', 'JAK': 'JAK1 JAK2',
        'p38': 'MAPK14 MAPKAPK2', 'JNK': 'MAPK8 MAPK9', 'PKC': 'PRKCA PRKCB PRKCD PRKCE', 'GSK3': 'GSK3A GSK3B', 'CDK4_6': 'CDK4 CDK6',
        'CDK': 'CDK1 CDK2 CDK7 CDK9 CDK12', 'DDR': 'ATM ATR PRKDC', 'CK2': 'CSNK2A1 CSNK2A2', 'ROCK': 'ROCK1 ROCK2', 'PAK': 'PAK1 PAK4',
        'CAMK2': 'CAMK2A CAMK2D', 'TAK1': 'MAP3K7', 'PKA': 'PRKACA', 'IKK': 'IKBKB', 'AURK_PLK': 'AURKA AURKB PLK1', 'RSK': 'RPS6KA1 RPS6KA3',
        'GF_RTK': 'MAPK1 MAPK3 RPS6KA1 AKT1', 'TNF': 'MAP3K7 IKBKB MAPK14 MAPKAPK2', 'DNA_damage': 'ATM ATR PRKDC', 'GPCR': 'MAPK1 MAPK3 PRKCA PRKACA',
        'PMA': 'PRKCA PRKCB PRKCD PRKCE', 'cytokine_other': 'JAK1 JAK2'}
TYR = {'EGFR', 'SRC', 'ABL1', 'FGFR1', 'FGFR2', 'MET', 'ALK', 'JAK1', 'JAK2'}
MP = pd.read_parquet(CA + 'motif_pct.parquet')
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'residue', 'psp_kinases_human'], low_memory=False).set_index('site')
BINS = [(0.5, 1), (1, 2), (2, 99)]


def build(kind):
    if kind == 'inh':
        X, E = load_exp('inh'); A = mask_binders(X, E); idx = X.index
    else:
        M = pd.read_csv(ACA + 'stim_meta.tsv', sep='\t', dtype={'study_pmid': str})
        drop = (M.control_raw == 'LPA') | (M.study_pmid == '29020638') | M.condition_raw.str.contains('p38i', na=False) | M.agent.isin(['cytarabine', 'methotrexate'])
        M = M[~drop]; Xc = pd.read_parquet(ACA + 'stim_cond.parquet')[M.conditionID]
        M['exp'] = M.study_pmid + '|' + M.cell_line_norm.astype(str) + '|' + M.cls
        X = Xc.T.groupby(M.exp.values).mean().T
        E = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), cls=('cls', 'first')).reindex(X.columns).reset_index()
        A = X.values.astype(np.float32); idx = X.index
    return A, E, idx


rows = []
prng = np.random.default_rng(7)
ALLK = [k for k in MP.columns if k not in TYR and MP[k].max() > 0][:]
for kind in ['inh', 'stim']:
    A, E, idx = build(kind)
    res = S.residue.reindex(idx).values
    psp = S.psp_kinases_human.reindex(idx).fillna('').values
    fin = ~np.isnan(A)
    for (c, kk0), draw in [(ck, dr) for ck in KMAP.items() for dr in range(5)]:
        c = c; kk = ' '.join(prng.choice([k for k in ALLK if k not in kk0.split()], len([k for k in kk0.split() if k not in TYR]) or 1, replace=False))
        ids = np.where(E.cls.values == c)[0]
        if len(ids) < 2 or E.study.iloc[ids].nunique() < 2: continue
        ks = kk.split(); kst = [k for k in ks if k not in TYR]; kty = [k for k in ks if k in TYR]
        mp = np.zeros(len(idx))
        sub = MP.reindex(idx)
        st_ = np.isin(res, ['S', 'T']); ty_ = res == 'Y'
        if kst: mp[st_] = sub.loc[st_, kst].max(axis=1).values
        if kty: mp[ty_] = sub.loc[ty_, kty].max(axis=1).values
        match = mp >= 90
        kset = set(ks); known = np.array([bool(set(p.split(';')) & kset) for p in psp])
        for i in ids:
            for j in ids:
                if E.study.iloc[i] == E.study.iloc[j]: continue
                both = fin[:, i] & fin[:, j]
                if both.sum() < 200: continue
                ai, aj = A[both, i], A[both, j]; mb = match[both]; kb = known[both]
                for sgn, lab in [(1, 'up'), (-1, 'down')]:
                    x, y = sgn * ai, sgn * aj; rep = y > 0.3; tot = len(y); rtot = rep.sum()
                    prim = (lab == 'down') if kind == 'inh' else (lab == 'up')
                    for lo, hi in [(1, 99)]:
                        d = (x > lo) & (x <= hi)
                        for split, sm in [('motif', mb)]:
                            for val in [True, False]:
                                dd = d & (sm == val); nd = dd.sum()
                                if nd < 5: continue
                                r = rep[dd].sum()
                                a_, b_, c_, d_ = r + .5, nd - r + .5, rtot - r + .5, (tot - nd) - (rtot - r) + .5
                                rows.append(dict(kind=kind, cls=c, draw=draw, kset=kk, si=E.study.iloc[i], sj=E.study.iloc[j], same_cell=E.cell.iloc[i] == E.cell.iloc[j],
                                                 prim='primary' if prim else 'reverse', bin=f'{lo}-{hi}', split=split, matched=val, n_disc=int(nd),
                                                 rep=r / nd, base=rtot / tot, lor=np.log(a_ * d_ / (b_ * c_))))
    print(kind, len(rows), flush=True)
R = pd.DataFrame(rows); R.to_csv(CA + 'direct_portability_PLACEBO_pairs.tsv.gz', sep='\t', index=False)
