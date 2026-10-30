"""H1: enrichment of axis sets among reverse vs primary responders, matched on residue x coverage x variance strata."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I25'); from common_i25 import *
OUT = '/workspace/results/b_I25/tables/'
SETS = [('CC_Mup', 'U_M'), ('CC_Mdown', 'U_M'), ('CC_MKI', 'U_MKI'), ('CC_union', 'U_CC'), ('CC_inter', 'U_INTER'),
        ('S_MET24', 'U_MET24'), ('S_MET24up', 'U_MET24'), ('S_MET24down', 'U_MET24'), ('S_MET30', 'U_MET30'), ('S_AMPKmotif', 'U_AMPKmotif'),
        ('S_CULT3D', 'U_CULT3D'), ('S_HIPPO', 'U_HIPPO'), ('S_SERUM', 'U_SERUM'), ('S_OIS', 'U_OIS')]
rows = []
for kind, incl in [('inh', False), ('inh', True), ('stim', False)]:
    A, E = load_test(kind, incl); AX = axes(A.index)
    strat = AX.stratum.fillna('NA').values
    prim = -1 if kind == 'inh' else 1
    for j, e in E.iterrows():
        x = A[e.exp].values; fin = ~np.isnan(x)
        for sname, uname in SETS:
            u = fin & AX[uname].values
            if u.sum() < 300: continue
            s = AX[sname].values[u]; xs = x[u]; st = strat[u]
            if s.sum() < 10: continue
            rate = pd.Series(s).groupby(st).mean()
            r_site = rate.reindex(st).values
            for lab, sg in [('primary', prim), ('reverse', -prim)]:
                resp = sg * xs > 1
                if resp.sum() < 10: continue
                rows.append(dict(kind=kind if not incl else 'inh_withCC', exp=e.exp, study=e.study, cell=e.cell, cls=e.cls, time_min=e.time_min, set=sname, dir=lab,
                                 n_resp=int(resp.sum()), O=float(s[resp].sum()), E=float(r_site[resp].sum()), O_raw_rate=s[resp].mean(), bg_rate=s.mean()))
R = pd.DataFrame(rows); R.to_csv(OUT + 'T1_H1_per_exp.tsv', sep='\t', index=False)

# aggregation in 02b_h1_agg.py
