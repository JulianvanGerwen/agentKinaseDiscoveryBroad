"""I10 step 1: experiment-level compendium matrix + perturbation group labels.
Experiment = study x cell line x treatment (doses of a dose series collapsed: mean over doses >= median dose of the series,
inhibitor dose series restricted to 300 nM - 10 uM as in I19). Time points of non-dose time series are kept separate for the
module analysis (exp_mod) but are averaged within the footprint experiments (exp_fp).
Each experiment vector is robust-z normalised (median-centred, divided by 1.4826*MAD) so studies are on one scale.
Outputs: cache/b_I10/{fp_meta.tsv, fp_sites.parquet, mod_meta.tsv, mod_sites.parquet}"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
sys.path.insert(0, "/workspace/scripts/b_I19")
OUT = "/workspace/cache/b_I10/"
# I19 inhibitor class map (agent -> class)
exec(open('/workspace/scripts/b_I19/02_build.py').read().split("A2C =")[0].split('OUT =')[1].split('\n', 1)[1])
A2C = {a: c for c, s in CLASS.items() for a in s.split()}
AXIS = {'MEK': 'MAPKaxis', 'ERK': 'MAPKaxis', 'RAF': 'MAPKaxis', 'PI3K': 'PI3Kaxis', 'PI3K_mTOR': 'PI3Kaxis', 'mTORcat': 'PI3Kaxis',
        'AKT': 'PI3Kaxis', 'PDPK1': 'PI3Kaxis', 'rapalog': 'PI3Kaxis', 'S6K': 'PI3Kaxis', 'CDK4_6': 'CDKaxis', 'CDK': 'CDKaxis'}

C = load.conditions()
C = C[~C.qc_flag.isin(['scale_suspect', 'duplicate_vector'])]
C = C[~C.perturbation_class.isin(['time_course'])]
C['sign'] = np.where(C.contrast_reversed.astype(bool), -1.0, 1.0)

def grp(r):
    pc, sc, st = r.perturbation_class, r.perturbation_subclass, r.stimulus
    if isinstance(r.background_agents, str) or r.is_combination: return None
    if pc == 'kinase_inhibitor':
        if sc not in ('inhibitor', 'inhibitor_reversed_contrast') or r.target_confidence not in ('high', 'medium'): return None
        return 'KI_' + A2C[r.kinase_inhibitors] if r.kinase_inhibitors in A2C else None
    if st in ('EGF', 'HB-EGF'): return 'EGF'
    if st == 'FGF': return 'FGF'
    if st == 'insulin': return 'insulin'
    if st == 'serum': return 'serum'
    if sc == 'starvation': return 'starvation'
    if sc == 'TNF': return 'TNF'
    if sc in ('interferon', 'IFNG_LPS', 'interleukin'): return 'JAKSTAT_cytokine'
    if sc in ('DNA_damage', 'UV'): return 'DNA_damage'
    if sc == 'replication_stress': return 'replication_stress'
    if sc == 'proteotoxic': return 'proteotoxic'
    if sc == 'oxidative': return 'oxidative'
    if sc == 'mito_depolarisation': return 'mito_depol'
    if sc == 'cell_cycle': return 'mitosis_cellcycle'
    if sc == 'GPCR_ligand': return 'GPCR_ligand'
    if sc == 'phosphatase_inhibitor': return 'phosphatase_inh'
    if sc == 'infection_or_coculture': return 'infection'
    if sc == 'immune_receptor': return 'immune_receptor'
    if sc == 'metabolic_drug': return 'metformin'
    return None
C['group'] = [grp(r) for r in C.itertuples()]
C['axis'] = C.group.map(lambda g: AXIS.get(g[3:]) if isinstance(g, str) and g.startswith('KI_') else None)
# inhibitor dose series in decryptM: restrict 300 nM..10 uM
isdec = (C.study_pmid == '36926954') & C.dose_series_id.notna()
keep_dose = ~isdec | C.dose_nM.between(300, 10000) | (C.perturbation_class != 'kinase_inhibitor')
# other dose series (non-inhibitor): keep doses >= series median
med = C.groupby('dose_series_id').dose_nM.transform('median')
keep_dose &= ~(C.dose_series_id.notna() & (C.perturbation_class != 'kinase_inhibitor')) | (C.dose_nM >= med)
C = C[keep_dose].copy()
C['treat'] = C.kinase_inhibitors.fillna(C.stimulus).fillna(C.agents).fillna(C.condition_raw.astype(str))
C['fp_exp'] = C.study_pmid + '|' + C.cell_line_norm.astype(str) + '|' + C.treat.astype(str) + '|' + C.background_agents.fillna('')
tm = C.time_min.fillna(-1).astype(int).astype(str)
C['mod_exp'] = np.where(C.dose_series_id.notna(), C.fp_exp + '|t' + tm, C.conditionID)
print(len(C), 'conditions;', C.fp_exp.nunique(), 'fp experiments;', C.mod_exp.nunique(), 'module experiments')

X = load.sites_matrix(C.conditionID.tolist())
X = X[X.notna().sum(1) >= 3]
X = X[C.conditionID] * C.sign.values

def agg_norm(key):
    W = X.T.groupby(C[key].values).mean().T
    m = W.median(); s = (W.sub(m)).abs().median() * 1.4826
    return (W.sub(m) / s).astype('float32')
for key, nm in [('fp_exp', 'fp'), ('mod_exp', 'mod')]:
    W = agg_norm(key)
    meta = C.groupby(key).agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), treat=('treat', 'first'), group=('group', 'first'),
                             paxis=('axis', 'first'), pclass=('perturbation_class', 'first'), subclass=('perturbation_subclass', 'first'),
                             n_cond=('conditionID', 'size'), time_min=('time_min', 'median'), dose_nM=('dose_nM', 'median'))
    meta['n_sites'] = W.notna().sum().reindex(meta.index)
    meta.index.name = 'exp'; meta.to_csv(OUT + f'{nm}_meta.tsv', sep='\t')
    W.to_parquet(OUT + f'{nm}_sites.parquet')
    print(nm, W.shape)
M = pd.read_csv(OUT + 'fp_meta.tsv', sep='\t', dtype={'study': str})
g = M.dropna(subset=['group']).groupby('group').agg(n_exp=('exp', 'size'), n_study=('study', 'nunique'), n_cell=('cell', 'nunique'))
g['eligible'] = (g.n_study >= 2) | (g.n_cell >= 3)
g.to_csv('/workspace/results/b_I10/tables/T0_groups.tsv', sep='\t'); print(g.sort_values('n_exp', ascending=False).to_string())
