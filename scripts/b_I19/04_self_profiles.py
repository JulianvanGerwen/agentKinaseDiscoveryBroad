"""Task 1b: robust self-hyperphosphorylation hits + feedback-vs-conformational dissection.
Robust hit (primary or binder tier): q<0.05 & Z>0, mean_lfc>0.4 vs DMSO, >=2 experiments, frac_exp_up>=0.6.
Also always-profiled controls: AKT1 S473/T308, BRAF/RAF1, MAP2K1/2 S218/S222, MAPK1/3 A-loop, MAPK14 T180/Y182.
For each site: values under own-binding drugs split by binding mode; under non-binding drugs of the same axis ('pathway-mates'); under unrelated drugs.
Classification (judgement): 'binding-associated' if up under binders but not under non-binding pathway-mates or unrelated drugs (|median|<0.25);
'feedback/network' if also up (median>0.4) under non-binding drugs that hit the same axis."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/b_I19/"; CA = OUT + 'cache/'
E = pd.read_csv(CA + 'experiments.tsv', sep='\t', dtype={'study': str}).fillna({'binders': '', 'primary': ''})
E['bset'] = E.binders.str.split(';').apply(lambda l: set(x for x in l if x))
X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]
R = pd.read_csv(OUT + 'tables/T1_self_sites_all.tsv', sep='\t')
S = load.sites().set_index('site')
# upstream axis for 'pathway-mate' test: classes whose inhibition should reduce input to / sit in same axis as kinase T (judgement)
UP = {'MAP2K1': ['ERK', 'RAF', 'EGFR', 'MET_ALK', 'RTKmulti'], 'MAP2K2': ['ERK', 'RAF', 'EGFR', 'MET_ALK', 'RTKmulti'],
      'MAPK1': ['MEK', 'RAF', 'EGFR'], 'MAPK3': ['MEK', 'RAF', 'EGFR'], 'BRAF': ['MEK', 'ERK'], 'RAF1': ['MEK', 'ERK'],
      'AKT1': ['PI3K', 'PI3K_mTOR', 'mTORcat', 'rapalog', 'PDPK1', 'S6K'], 'AKT2': ['PI3K', 'PI3K_mTOR', 'mTORcat', 'PDPK1'],
      'EGFR': ['MEK', 'ERK', 'RAF', 'PI3K', 'AKT'], 'ERBB2': ['MEK', 'ERK', 'PI3K', 'AKT'], 'MAPK14': ['TAK1'], 'CDK1': ['CDK4_6'],
      'PRKCA': ['PDPK1', 'mTORcat'], 'SRC': ['EGFR', 'MET_ALK'], 'LCK': [], 'PAK1': [], 'PAK4': []}
hits = R[(R.q < 0.05) & (R.Z > 0) & (R.mean_lfc > 0.4) & (R.n_exp >= 2) & (R.frac_exp_up >= 0.6)].copy()
ctrl_sites = {'AKT1': ['S473', 'T308'], 'MAP2K1': ['S218', 'S222'], 'MAP2K2': ['S222', 'S226'], 'MAPK1': ['T185', 'Y187'], 'MAPK3': ['T202', 'Y204'],
              'BRAF': ['S365', 'S446', 'T599', 'S602', 'S729'], 'RAF1': ['S259', 'S338', 'S621', 'S642'], 'MAPK14': ['T180', 'Y182'], 'EGFR': ['Y1197', 'Y1092', 'Y1172', 'S1166', 'T693'], 'AKT2': ['S474']}
ctrl = []
for g, ss in ctrl_sites.items():
    for s in S.index[(S.gene == g) & (S.residue + S.position.astype(str)).isin(ss)]:
        if s in X.index: ctrl.append((s, g))
todo = list(dict.fromkeys([(r.site, r.gene) for r in hits.itertuples()] + ctrl))
rows = []
for s, g in todo:
    x = X.loc[s]
    own = E[E.bset.apply(lambda b: g in b)]
    mates = E[E.cls.isin(UP.get(g, [])) & ~E.bset.apply(lambda b: g in b)]
    unrel = E[~E.cls.isin(UP.get(g, [])) & ~E.bset.apply(lambda b: g in b)]
    def summ(df):
        v = x[df.exp].dropna(); return (len(v), float(v.median()) if len(v) else np.nan, float((v > 0.4).mean()) if len(v) else np.nan)
    d = dict(site=s, gene=g, residue_pos=s.split('_')[1])
    for nm, df in [('own', own), ('own_ATPcomp', own[~own['mode'].isin(['allosteric', 'non-ATP'])]), ('own_allosteric', own[own['mode'].isin(['allosteric', 'non-ATP'])]), ('mates', mates), ('unrelated', unrel)]:
        n, med, fu = summ(df); d.update({f'{nm}_n': n, f'{nm}_med': med, f'{nm}_fracup': fu})
    d['own_detail'] = '; '.join(f"{r.agent}/{r.cell_line}:{x[r.exp]:.2f}" for r in own.itertuples() if pd.notna(x[r.exp]))
    d['mates_detail'] = '; '.join(f"{r.agent}/{r.cell_line}:{x[r.exp]:.2f}" for r in mates.itertuples() if pd.notna(x[r.exp]))
    om = d['own_med']; mm = d['mates_med']; um = d['unrelated_med']
    if om > 0.4 and (pd.notna(mm) and mm > 0.4): cl = 'feedback/network (also up under non-binding axis drugs)'
    elif om > 0.4 and (pd.isna(mm) or abs(mm) < 0.25) and abs(um) < 0.25: cl = 'binding-associated (own-drug specific)'
    elif om > 0.4: cl = 'up under own drug; ambiguous'
    else: cl = 'not up'
    d['call'] = cl; rows.append(d)
P = pd.DataFrame(rows).merge(R[R.tier == 'primary'][['site', 'aloop', 'psp_function', 'Z', 'q']], on='site', how='left')
P.to_csv(OUT + 'tables/T1_self_hits_profiles.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_colwidth', 200); pd.set_option('display.max_rows', 200)
print(P[['site', 'gene', 'residue_pos', 'own_n', 'own_med', 'own_ATPcomp_med', 'own_allosteric_n', 'own_allosteric_med', 'mates_n', 'mates_med', 'unrelated_n', 'unrelated_med', 'unrelated_fracup', 'call']].to_string())
for r in P.itertuples():
    print(r.site, r.gene, '| OWN:', r.own_detail[:400], '| MATES:', str(r.mates_detail)[:300])
