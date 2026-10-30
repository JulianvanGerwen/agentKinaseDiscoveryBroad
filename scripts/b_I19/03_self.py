"""Task 1: self-hyperphosphorylation catalogue.
For each kinase T that is a primary target (tier 'primary') or Klaeger/Hijazi binder (tier 'binder') of a drug, and each measured own phosphosite,
compare the site's log2FC under T-binding drugs with (i) zero (vs DMSO) and (ii) the same study's experiments whose drugs do NOT bind T.
Per experiment z = (x - median_others_in_study) / sd_others (robust, floored); Stouffer across experiments within study -> study z; Stouffer across studies.
Feedback vs conformational: also compute the site's response under drugs of the same pathway that do not bind T ('pathway-mate' drugs)."""
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/b_I19/"; CA = OUT + 'cache/'
E = pd.read_csv(CA + 'experiments.tsv', sep='\t', dtype={'study': str}).fillna({'binders': '', 'primary': '', 'klaeger1000': '', 'hijazi': ''})
X = pd.read_parquet(CA + 'site_exp.parquet')[E.exp]
S = load.sites().set_index('site')
al = pd.read_csv(CA + 'aloop.tsv', sep='\t').drop_duplicates('uniprot').set_index('uniprot')
reg = pd.read_csv('/workspace/data/Regulatory_sites', sep='\t', skiprows=3, low_memory=False)
reg = reg[(reg.ORGANISM == 'human') & reg.MOD_RSD.str.endswith('-p', na=False)]
reg['site'] = reg.ACC_ID + '_' + reg.MOD_RSD.str[:-2]
regf = reg.groupby('site').ON_FUNCTION.first()
# pathway membership (for 'pathway-mate' feedback test): classes considered to be in the same linear axis
AXIS = {'MEK': 'MAPK', 'ERK': 'MAPK', 'RAF': 'MAPK', 'EGFR': 'RTK', 'MET_ALK': 'RTK', 'RTKmulti': 'RTK', 'FGFR': 'RTK',
        'PI3K': 'PI3K', 'PI3K_mTOR': 'PI3K', 'mTORcat': 'PI3K', 'rapalog': 'PI3K', 'AKT': 'PI3K', 'PDPK1': 'PI3K', 'S6K': 'PI3K'}
E['bset'] = E.binders.str.split(';').apply(lambda l: set(x for x in l if x))
E['pset'] = E.primary.str.split(';').apply(lambda l: set(x for x in l if x))
gene = S.gene.reindex(X.index)
targets = sorted(set().union(*E.bset))
X = X[gene.isin(targets).values]; gene = gene.reindex(X.index)
print('own-kinase sites measured:', len(X))
st_of = E.set_index('exp').study
rows = []
for s, g in gene.items():
    x = X.loc[s]
    for tier, col in [('primary', 'pset'), ('binder', 'bset')]:
        own = E[E[col].apply(lambda b: g in b) & x.reindex(E.exp).notna().values]
        if tier == 'binder': own = own[~own.pset.apply(lambda b: g in b)]
        if len(own) == 0: continue
        zs, recs = [], []
        for st, ow in own.groupby('study'):
            oth = E[(E.study == st) & ~E.bset.apply(lambda b: g in b)]
            xo = x[oth.exp].dropna()
            xv = x[ow.exp].values
            if len(xo) >= 3:
                med = xo.median(); sd = max(1.4826 * (xo - med).abs().median(), 0.15)
                ze = (xv - med) / sd
                zs.append(ze.sum() / np.sqrt(len(ze)))
                d = xv.mean() - med
            else:
                d = np.nan
            recs.append((st, len(ow), xv.mean(), d))
        n_st = len(recs)
        zst = np.array(zs)
        Z = zst.sum() / np.sqrt(len(zst)) if len(zst) else np.nan
        u = S.loc[s, 'uniprot'] if s in S.index else None; pos = S.loc[s, 'position']
        aloop = (u in al.index) and (al.loc[u, 'dfg_pos'] <= pos <= al.loc[u, 'ape_pos'])
        rows.append(dict(site=s, gene=g, tier=tier, residue=S.loc[s, 'residue'], position=pos, aloop=aloop, psp_function=regf.get(s, np.nan),
                         n_exp=len(own), n_studies=n_st, n_agents=own.agent.nunique(), n_cell=own.cell_line.nunique(),
                         agents=';'.join(sorted(own.agent.unique())), modes=';'.join(sorted(own['mode'].unique())),
                         mean_lfc=float(np.nanmean([r[2] for r in recs])), frac_exp_up=float((x[own.exp] > 0).mean()),
                         mean_delta_vs_other=float(np.nanmean([r[3] for r in recs])) if any(~np.isnan([r[3] for r in recs])) else np.nan,
                         n_studies_with_ref=len(zst), n_studies_up=int((zst > 0).sum()), Z=Z,
                         per_study=';'.join(f'{r[0]}:{r[1]}:{r[2]:.2f}' for r in recs)))
R = pd.DataFrame(rows)
R['p'] = 2 * stats.norm.sf(R.Z.abs())
for tier in ['primary', 'binder']:
    m = R.tier == tier
    R.loc[m & R.p.notna(), 'q'] = stats.false_discovery_control(R.loc[m & R.p.notna(), 'p'])
R.to_csv(OUT + 'tables/T1_self_sites_all.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200); pd.set_option('display.max_colwidth', 50)
cols = ['site', 'gene', 'residue', 'position', 'aloop', 'psp_function', 'n_exp', 'n_studies', 'n_agents', 'n_cell', 'mean_lfc', 'mean_delta_vs_other', 'frac_exp_up', 'n_studies_up', 'n_studies_with_ref', 'Z', 'q', 'agents']
print(R[(R.tier == 'primary')].sort_values('Z', ascending=False).head(40)[cols].to_string())
print(R[(R.tier == 'binder')].sort_values('Z', ascending=False).head(40)[cols].to_string())
