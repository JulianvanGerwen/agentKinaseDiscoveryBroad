"""Cell-line contrast: under growth-factor stimulation and RTK (EGFR/ERBB2/MET/FGFR/...) inhibition, do RTK pY, adaptor pY,
ERK activation-loop phosphorylation co-move (transitive)? Same IV decomposition across conditions."""
from common_a import zs, meta_re, OUT
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, '/workspace/scripts/infra_compendium'); import load
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 40)
C = load.conditions(); S = load.sites()
RTKG = ['EGFR', 'ERBB2', 'ERBB3', 'ERBB4', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'PDGFRA', 'PDGFRB', 'IGF1R', 'INSR', 'KIT', 'AXL',
        'NTRK1', 'NTRK2', 'RET', 'ALK', 'ROS1', 'EPHA2', 'DDR1', 'FLT3', 'KDR', 'MST1R']
RTK_ACT = {'EGFR', 'ERBB2', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'PDGFRA', 'PDGFRB', 'IGF1R', 'INSR', 'ALK', 'RET', 'KIT', 'FLT3', 'AXL', 'ERBB3'}
ADAPT = ['P29353_Y349', 'P29353_Y350', 'P29353_Y427', 'Q13480_Y627', 'Q13480_Y659', 'Q06124_Y542', 'Q06124_Y580']
ERK = ['P27361_T202', 'P27361_Y204', 'P28482_T185', 'P28482_Y187']
ERK_I24 = ERK
rtk_sites = list(S[S.gene.isin(RTKG) & (S.residue == 'Y')].site)
# conditions: growth factor stimulation + RTK inhibitor (single agent, primary target in RTK set)
def is_rtki(t):
    if not isinstance(t, str): return False
    return any(g in RTK_ACT for g in t.replace(';', ',').replace('|', ',').split(','))
gf = C[(C.perturbation_class == 'growth_factor') & ~C.contrast_reversed.fillna(False).astype(bool)]
ki = C[(C.perturbation_class == 'kinase_inhibitor') & C.inhibitor_primary_targets.apply(is_rtki) & (C.is_combination != True)]
sel = pd.concat([gf.assign(kind='GF'), ki.assign(kind='RTKi')])
sel = sel[sel.qc_flag.isna() | (sel.qc_flag == 'ok')]
X = load.sites_matrix(sel.conditionID, sites=rtk_sites + ADAPT + ERK)
def comp(sites):
    ss = [s for s in sites if s in X.index]; return X.loc[ss].mean(0, skipna=True), X.loc[ss].notna().sum(0)
D = sel.set_index('conditionID')[['study_pmid', 'kind', 'stimulus', 'agents', 'inhibitor_primary_targets', 'cell_line_norm', 'time_min', 'dose_nM']].copy()
for k, s in [('RTK', rtk_sites), ('ADAPT', ADAPT), ('ERK', ERK)]:
    v, n = comp(s); D[k] = v.reindex(D.index); D[k + '_n'] = n.reindex(D.index)
D = D.dropna(subset=['RTK', 'ADAPT', 'ERK'])
D.to_csv(f'{OUT}/tables/A4_compendium_conditions.tsv', sep='\t')
print(D.groupby('kind').size(), D.study_pmid.nunique(), 'studies')
print(D.groupby(['kind', 'study_pmid']).size().to_string())
# sign concordance
for k, g in D.groupby('kind'):
    print(k, 'n', len(g), 'median log2FC RTK %.2f ADAPT %.2f ERK %.2f' % (g.RTK.median(), g.ADAPT.median(), g.ERK.median()),
          '| frac same sign RTK&ADAPT %.2f, ADAPT&ERK %.2f, RTK&ERK %.2f' % (np.mean(np.sign(g.RTK) == np.sign(g.ADAPT)),
          np.mean(np.sign(g.ADAPT) == np.sign(g.ERK)), np.mean(np.sign(g.RTK) == np.sign(g.ERK))))
# correlations across all conditions (pooled, Spearman) and within-study-demeaned
def gains(d, lab):
    out = {}
    for x, y in [('RTK', 'ADAPT'), ('ADAPT', 'ERK'), ('RTK', 'ERK')]:
        out[f'{x}>{y}'] = stats.spearmanr(d[x], d[y])[0]
    z = d[['RTK', 'ADAPT', 'ERK']].apply(zs)
    fs = np.polyfit(z.RTK, z.ADAPT, 1)[0]; rf = np.polyfit(z.RTK, z.ERK, 1)[0]
    Xm = np.c_[np.ones(len(z)), z.ADAPT, z.RTK]; oth = np.linalg.lstsq(Xm, z.ERK, rcond=None)[0][1]
    rng = np.random.default_rng(0); B = []
    for _ in range(2000):
        i = rng.integers(0, len(z), len(z)); zz = z.iloc[i]
        f_ = np.polyfit(zz.RTK, zz.ADAPT, 1)[0]; r_ = np.polyfit(zz.RTK, zz.ERK, 1)[0]
        Xb = np.c_[np.ones(len(zz)), zz.ADAPT, zz.RTK]; o_ = np.linalg.lstsq(Xb, zz.ERK.values, rcond=None)[0][1]
        B.append([r_ / f_, o_, r_ / f_ - o_])
    B = np.array(B)
    out.update(lab=lab, n=len(d), studies=d.study_pmid.nunique(), iv=rf / fs, iv_lo=np.percentile(B[:, 0], 2.5), iv_hi=np.percentile(B[:, 0], 97.5),
               other=oth, diff=rf / fs - oth, diff_lo=np.percentile(B[:, 2], 2.5), diff_hi=np.percentile(B[:, 2], 97.5), product=out['RTK>ADAPT'] * out['ADAPT>ERK'])
    return out
rows = [gains(D, 'all'), gains(D[D.kind == 'GF'], 'GF'), gains(D[D.kind == 'RTKi'], 'RTKi')]
W = D.copy()
for k in ['RTK', 'ADAPT', 'ERK']: W[k] = W[k] - W.groupby('study_pmid')[k].transform('mean')
W = W[W.groupby('study_pmid').RTK.transform('size') >= 3]
rows.append(gains(W, 'all_within_study'))
G = pd.DataFrame(rows); G.to_csv(f'{OUT}/tables/A4_compendium_transitivity.tsv', sep='\t', index=False); print(G.round(3).to_string())

# refined: EGFR/ERBB2 axis only. RTK pY = EGFR/ERBB2 pY sites; conditions = EGF/HB-EGF/serum stimulation or EGFR/ERBB2-targeted inhibitor, dose>=100 nM
eg_sites = list(S[S.gene.isin(['EGFR', 'ERBB2']) & (S.residue == 'Y')].site)
def egfr_t(t): return isinstance(t, str) and any(g in ('EGFR', 'ERBB2') for g in t.split(';'))
sel2 = pd.concat([gf[gf.stimulus.isin(['EGF', 'HB-EGF', 'serum'])].assign(kind='EGF'),
                  ki[ki.inhibitor_primary_targets.apply(egfr_t) & ((ki.dose_nM >= 100) | ki.dose_nM.isna())].assign(kind='EGFRi')])
sel2 = sel2[sel2.qc_flag.isna() | (sel2.qc_flag == 'ok')]
X = load.sites_matrix(sel2.conditionID, sites=eg_sites + ADAPT + ERK)
D2 = sel2.set_index('conditionID')[['study_pmid', 'kind', 'stimulus', 'agents', 'cell_line_norm', 'time_min', 'dose_nM']].copy()
for k, s in [('RTK', eg_sites), ('ADAPT', ADAPT), ('ERK', ERK)]:
    v, n = comp(s); D2[k] = v.reindex(D2.index); D2[k + '_n'] = n.reindex(D2.index)
D2 = D2.dropna(subset=['RTK', 'ADAPT', 'ERK']); D2.to_csv(f'{OUT}/tables/A4b_compendium_EGFR_conditions.tsv', sep='\t')
print(D2.groupby(['kind', 'study_pmid']).size().to_string())
for k, g in D2.groupby('kind'):
    print(k, 'n', len(g), 'median log2FC RTK %.2f ADAPT %.2f ERK %.2f' % (g.RTK.median(), g.ADAPT.median(), g.ERK.median()),
          '| same sign RTK&ADAPT %.2f, ADAPT&ERK %.2f, RTK&ERK %.2f' % (np.mean(np.sign(g.RTK) == np.sign(g.ADAPT)),
          np.mean(np.sign(g.ADAPT) == np.sign(g.ERK)), np.mean(np.sign(g.RTK) == np.sign(g.ERK))))
rows = [gains(D2, 'EGFRaxis_all'), gains(D2[D2.kind == 'EGFRi'], 'EGFRaxis_EGFRi')]
if (D2.kind == 'EGF').sum() >= 8: rows.append(gains(D2[D2.kind == 'EGF'], 'EGFRaxis_EGF'))
W = D2.copy()
for k in ['RTK', 'ADAPT', 'ERK']: W[k] = W[k] - W.groupby('study_pmid')[k].transform('mean')
W = W[W.groupby('study_pmid').RTK.transform('size') >= 3]; rows.append(gains(W, 'EGFRaxis_within_study'))
G2 = pd.DataFrame(rows); G2.to_csv(f'{OUT}/tables/A4b_compendium_EGFR_transitivity.tsv', sep='\t', index=False); print(G2.round(3).to_string())
