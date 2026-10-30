"""Audit I19 - independent build.
Inhibitor conditions: same inclusion rules as stated by I19 (re-implemented, not imported), kept at condition level
(for time/dose matching) and aggregated to experiments (study x cell x agent).
Stimulation conditions (mirror test): RTK ligands / serum, TNF, DNA-damage stress, GPCR ligands, phosphatase inhibitors.
Outputs: results/audit_I19/cache/{inh_cond,stim_cond}.parquet (sites x conditions, sign-corrected), meta tsvs."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/audit_I19/cache/"
# agent -> class (copied judgement table from b_I19/02_build.py so that experiment sets are comparable)
CLASS = {'MEK': 'PD0325901 selumetinib refametinib U0126 trametinib cobimetinib PD184352', 'ERK': 'ravoxertinib ulixertinib ERKi_generic ERKi_SCH772984',
 'RAF': 'vemurafenib dabrafenib SB590885 LY3009120 sorafenib', 'EGFR': 'lapatinib gefitinib erlotinib afatinib CUDC101 PD153035 AZD3759 AG1478 PD168393',
 'PI3K': 'pictilisib HS173 PIK294 idelalisib TGX221 wortmannin LY294002 PI3Ka_inh2 AZD6482', 'PI3K_mTOR': 'dactolisib PI103 omipalisib',
 'mTORcat': 'AZD8055 KU0063794 torin1 torkinib', 'rapalog': 'rapamycin temsirolimus everolimus',
 'AKT': 'MK2206 GSK690693 capivasertib ipatasertib akti_akt12 aktVIII AT13148', 'PDPK1': 'GSK2334470 BX912', 'S6K': 'PF4708671 LY2584702 DG2',
 'SRC_ABL': 'dasatinib bosutinib imatinib', 'FGFR': 'AZD4547', 'MET_ALK': 'crizotinib PF04217903 JNJ38877605 foretinib lorlatinib ceritinib TAE684',
 'JAK': 'tofacitinib AZD1480', 'p38': 'PH797804 TAK715 SB203580 BIRB796 SB202190 SB239063', 'JNK': 'SP600125 JNKIN8 tanzisertib',
 'PKC': 'BIM1 Go6983 Go6976 enzastaurin', 'GSK3': 'tideglusib CHIR99021 LY2090314 GSK3i_BIO SB216763', 'CDK4_6': 'ribociclib abemaciclib',
 'CDK': 'AZD5438 RO3306 flavopiridol THZ531 CDK12IN3 kenpaullone cortistatinA', 'DDR': 'VE821 AZ20 ceralasertib NU7026 NU7441 KU60019', 'CK2': 'CX4945 TBCA',
 'ROCK': 'ripasudil belumosudil H1152 Y27632', 'PAK': 'PF3758309 FRAX486', 'CAMK2': 'KN62 KN93', 'TAK1': '7Oxo',
 'RTKmulti': 'nintedanib amuvatinib quizartinib pexidartinib linsitinib', 'PKA': 'H89', 'IKK': 'SC514', 'AURK_PLK': 'MLN8054 AZD1152 VX680 BI2536 BI4834',
 'RSK': 'BID1870', 'other': 'LRRK2IN1 ibrutinib'}
A2C = {a: c for c, s in CLASS.items() for a in s.split()}
C = load.conditions(); D = load.drug_targets().set_index('agent')
bad_qc = C.qc_flag.isin(['scale_suspect', 'duplicate_vector'])
ki = C[(C.perturbation_class == 'kinase_inhibitor') & C.target_confidence.isin(['high', 'medium']) & ~bad_qc
       & C.perturbation_subclass.isin(['inhibitor', 'inhibitor_reversed_contrast']) & C.background_agents.isna()
       & ~C.is_combination.astype(bool) & C.kinase_inhibitors.isin(A2C)].copy()
ki = ki[~(ki.dose_series_id.notna() & (ki.study_pmid == '36926954')) | ki.dose_nM.between(300, 10000)]
ki['cls'] = ki.kinase_inhibitors.map(A2C); ki['agent'] = ki.kinase_inhibitors
def binders(a):
    if a not in D.index: return ''
    r = D.loc[a]; return ';'.join(str(r[c]) for c in ['inhibitor_primary_targets' if 'inhibitor_primary_targets' in D.columns else 'primary_targets', 'klaeger_kd_le1000nM', 'hijazi_targets'] if c in D.columns and pd.notna(r[c]))
ki['binders'] = ki.agent.map(binders)
# stimulations
DNA = {'ionising radiation', 'UV'}; DNAd = {'etoposide', 'hydroxyurea', 'neocarzinostatin', 'doxorubicin', 'bleomycin', 'cisplatin', 'carboplatin', '4NQO', 'temozolomide', 'cytarabine', 'methotrexate'}
GPCR = {'endothelin', 'LPA', 'thrombin', 'PGE2', 'angiotensin II', 'ADP'}
def stim_cls(r):
    if r.perturbation_class == 'growth_factor': return 'GF_RTK'
    if r.stimulus == 'TNF': return 'TNF'
    if r.perturbation_class == 'cytokine': return 'cytokine_other'
    if r.stimulus in DNA or (isinstance(r.agents, str) and r.agents in DNAd and r.perturbation_class == 'stress'): return 'DNA_damage'
    if r.stimulus in GPCR: return 'GPCR'
    if r.perturbation_subclass == 'phosphatase_inhibitor': return 'PPase_inh'
    if r.stimulus == 'PMA': return 'PMA'
    if isinstance(r.agents, str) and r.agents in {'bortezomib', 'carfilzomib'}: return 'proteasome'
    return None
st = C[~bad_qc & C.background_agents.isna() & ~C.is_combination.astype(bool)].copy()
st['cls'] = [stim_cls(r) for r in st.itertuples()]
st = st[st.cls.notna()].copy(); st['agent'] = st.stimulus.fillna(st.agents); st['binders'] = ''
for df in (ki, st):
    df['sign'] = np.where(df.contrast_reversed.astype(bool), -1.0, 1.0)
    df['exp'] = df.study_pmid + '|' + df.cell_line_norm.astype(str) + '|' + df.agent.astype(str)
cols = ['conditionID', 'study_pmid', 'cell_line_norm', 'agent', 'cls', 'exp', 'sign', 'dose_nM', 'time_min', 'binders', 'control_raw', 'condition_raw']
ki[cols].to_csv(OUT + 'inh_meta.tsv', sep='\t', index=False); st[cols].to_csv(OUT + 'stim_meta.tsv', sep='\t', index=False)
for name, df in [('inh', ki), ('stim', st)]:
    X = load.sites_matrix(df.conditionID.tolist())
    X = X * df.sign.values
    X = X[X.notna().any(axis=1)]
    X.astype('float32').to_parquet(OUT + f'{name}_cond.parquet')
    Xe = X.T.groupby(df.exp.values).mean().T
    Xe.astype('float32').to_parquet(OUT + f'{name}_exp.parquet')
    print(name, len(df), 'conds', df.exp.nunique(), 'exps', df.study_pmid.nunique(), 'studies', X.shape)
print(st.groupby('cls').agg(n=('conditionID', 'size'), nexp=('exp', 'nunique'), ns=('study_pmid', 'nunique')))
print(ki.groupby('cls').agg(n=('conditionID', 'size'), nexp=('exp', 'nunique'), ns=('study_pmid', 'nunique')).sort_values('ns'))
