"""Task 4: cross-cohort meta of reduced-form / Wald IV per kinase x site; discovery and unsupported-PSP lists.
Meta: IVW of Wald ratio beta_iv = b_rf / b_fs (SE = SE_rf/|b_fs|, strong-instrument approx) over cohorts where kinase F>10.
Discovery: meta q<0.05 (BH within kinase), beta>0, nominal p<0.05 & beta>0 in >=2 cohorts, motif pct>=90, not PSP for kinase,
site not on kinase's chromosome arm; robust if neighbour-adjusted t>1 in the majority of supporting cohorts.
Null: same pipeline with the motif filter swapped for an unrelated kinase's motif (different group) -> empirical FDR of motif step."""
import numpy as np, pandas as pd
from scipy import stats
from lib import *
rng = np.random.default_rng(1)
ann = pd.read_parquet(f'{CACHE}/gene_annot.parquet').drop_duplicates('uniprot').set_index('uniprot')
psp = pd.read_parquet(f'{CACHE}/psp_human.parquet')
psp['kin'] = psp.KIN_ACC_ID.str.split('-').str[0]
psp['site_key'] = psp.SUB_ACC_ID.str.split('-').str[0] + '_' + psp.SUB_MOD_RSD
psp_pairs = set(zip(psp.kin, psp.site_key)); psp_any = set(psp.site_key)
motif = pd.read_parquet(f'{CACHE}/motif_pct.parquet')
kl = pd.read_csv('<HACKATHON_ROOT>/kinaseSubMain/results/labels/kinase_list.tsv', sep='\t')
grp = kl.drop_duplicates('gene').set_index('gene').manning_group.to_dict()
sym = ann.symbol.to_dict()

long = []
for C in COHORTS:
    d = np.load(f'{CACHE}/{C}/rf.npz', allow_pickle=True)
    t, b, tnb, tfo, tnv = d['t_rf'], d['b_rf'], d['t_nb'], d['t_fo'], d['t_nv']
    se = np.abs(b / t)
    biv = b / d['bfs'][:, None]; seiv = se / np.abs(d['bfs'][:, None])
    K, S = t.shape
    kk = np.repeat(d['kin'], S); ss = np.tile(d['sites'], K)
    long.append(pd.DataFrame({'cohort': C, 'kin': kk, 'site': ss, 't': t.ravel(), 'biv': biv.ravel(), 'seiv': seiv.ravel(),
                              't_nb': tnb.ravel(), 't_fo': tfo.ravel(), 't_nv': tnv.ravel(), 'n': np.tile(d['n_s'], K)}))
L = pd.concat(long, ignore_index=True)
L = L[np.isfinite(L.t) & np.isfinite(L.seiv) & (L.seiv > 0)]
L['p'] = 2 * stats.t.sf(np.abs(L.t), L.n - 3)
L['w'] = 1 / L.seiv ** 2
L['bw'] = L.biv * L.w
L['pos_nom'] = (L.p < 0.05) & (L.t > 0)
L['neg_nom'] = (L.p < 0.05) & (L.t < 0)
L['nb_ok'] = L.t_nb > 1
L['nb_avail'] = np.isfinite(L.t_nb)
g = L.groupby(['kin', 'site'])
M = g.agg(k=('t', 'size'), sw=('w', 'sum'), sbw=('bw', 'sum'), n_pos=('pos_nom', 'sum'), n_neg=('neg_nom', 'sum'),
          n_nb_ok=('nb_ok', 'sum'), n_nb_av=('nb_avail', 'sum'), mean_t_nv=('t_nv', 'mean'), mean_t_fo=('t_fo', 'mean'))
cp = L[L.pos_nom].groupby(['kin', 'site']).cohort.agg(','.join).rename('cohorts_pos')
M = M.join(cp)
M['beta_iv'] = M.sbw / M.sw; M['se_iv'] = 1 / np.sqrt(M.sw); M['z'] = M.beta_iv / M.se_iv
M['p_meta'] = 2 * stats.norm.sf(np.abs(M.z))
M = M.reset_index()
M['q_meta'] = M.groupby('kin').p_meta.transform(lambda p: pd.Series(bh(p.values), index=p.index))
M['kgene'] = M.kin.map(sym)
M['site_u'] = M.site.str.split('_').str[0]
M['sgene'] = M.site_u.map(sym)
M['cis'] = (M.site_u == M.kin) | (M.site_u.map(ann.arm) == M.kin.map(ann.arm))
M['psp_edge'] = pd.MultiIndex.from_arrays([M.kin, M.site]).isin(list(psp_pairs))
M['psp_anysite'] = M.site.isin(psp_any)
mv = motif.stack().rename('pct').reset_index().rename(columns={'level_0': 'site', 'level_1': 'kgene'}) if False else None
mot = motif.reindex(columns=sorted(set(M.kgene.dropna()) & set(motif.columns)))
ms = mot.stack().rename('motif_pct'); ms.index.names = ['site', 'kgene']
M = M.merge(ms.reset_index(), on=['site', 'kgene'], how='left')
M.to_parquet(f'{RES}/tables/T3_meta_kinase_site.parquet')

base = (M.k >= 2) & (M.q_meta < 0.05) & (M.beta_iv > 0) & (M.n_pos >= 2) & ~M.cis
disc = M[base & (M.motif_pct >= 90) & ~M.psp_edge].copy()
disc['nb_robust'] = (disc.n_nb_av == 0) | (disc.n_nb_ok >= np.ceil(disc.n_nb_av / 2))
disc = disc.sort_values('p_meta')
disc.to_csv(f'{RES}/tables/T4_discovery_edges.tsv', sep='\t', index=False)
# empirical motif null: replace kinase motif with a random kinase of different Manning group
kg = M.kgene.dropna().unique(); allm = [c for c in motif.columns]
nullc = []
for r in range(20):
    swap = {}
    for k in kg:
        cands = [c for c in allm if grp.get(c, 'x') != grp.get(k, 'y')]
        swap[k] = rng.choice(cands)
    sub = M[base & ~M.psp_edge]
    pv = [motif.at[s, swap[k]] if (k in swap and s in motif.index) else np.nan for s, k in zip(sub.site, sub.kgene)]
    nullc.append((np.array(pv) >= 90).sum())
print('base hits (meta q<0.05, >=2 cohorts pos, non-cis, non-PSP):', (base & ~M.psp_edge).sum())
print('  of which own-motif>=90:', len(disc), ' null (other-group motif) mean', np.mean(nullc), 'sd', np.std(nullc))
print('  base hits with motif scored:', (base & ~M.psp_edge & M.motif_pct.notna()).sum())
print('  frac own motif>=90 among scored:', (disc.shape[0]) / max(1, (base & ~M.psp_edge & M.motif_pct.notna()).sum()))
print('PSP edges among base hits:', (base & M.psp_edge).sum(), 'of tested PSP edges', (M.psp_edge & (M.k >= 2) & ~M.cis).sum(),
      ' vs base rate', base.mean())
print('nb-robust discoveries:', disc.nb_robust.sum())
print(disc.head(40)[['kgene', 'site', 'sgene', 'k', 'n_pos', 'cohorts_pos', 'beta_iv', 'se_iv', 'z', 'q_meta', 'motif_pct', 'psp_anysite', 'n_nb_ok', 'n_nb_av', 'mean_t_nv', 'mean_t_fo']].to_string())
# unsupported PSP edges: well powered (k>=3), meta IV upper 95% CI < 0.2
uns = M[M.psp_edge & (M.k >= 3) & ~M.cis].copy()
uns['upper'] = uns.beta_iv + 1.96 * uns.se_iv
uns['lower'] = uns.beta_iv - 1.96 * uns.se_iv
uns['status'] = np.select([uns.lower > 0, uns.upper < 0, uns.upper < 0.2], ['supported', 'contradicted(neg)', 'powered_null'], 'uninformative')
uns.sort_values('z').to_csv(f'{RES}/tables/T5_psp_edges_status.tsv', sep='\t', index=False)
print(uns.status.value_counts())
print(uns[uns.status != 'uninformative'].sort_values('z')[['kgene', 'site', 'sgene', 'k', 'beta_iv', 'se_iv', 'z', 'n_pos', 'n_neg', 'status']].to_string())
