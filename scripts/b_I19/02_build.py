"""Build experiment-level tables for I19.
Experiment = study x cell line x agent (x time kept as mean over time points/doses).
Inhibitor conditions: kinase_inhibitor class, target_confidence high/medium, qc ok/low_coverage? (exclude scale_suspect/duplicate),
exclude inhibitor_vs_other_inhibitor, inhibitor+stimulus, background agents, combinations; flip contrast_reversed.
decryptM/dose series: keep doses 300 nM..10 uM (avoid sub-target and >10uM toxic doses).
Outputs in results/b_I19/cache/."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/b_I19/cache/"
CLASS = {  # agent -> pathway class (judgement; based on curated primary targets)
 'MEK': 'PD0325901 selumetinib refametinib U0126 trametinib cobimetinib PD184352',
 'ERK': 'ravoxertinib ulixertinib ERKi_generic ERKi_SCH772984',
 'RAF': 'vemurafenib dabrafenib SB590885 LY3009120 sorafenib',
 'EGFR': 'lapatinib gefitinib erlotinib afatinib CUDC101 PD153035 AZD3759 AG1478 PD168393',
 'PI3K': 'pictilisib HS173 PIK294 idelalisib TGX221 wortmannin LY294002 PI3Ka_inh2 AZD6482',
 'PI3K_mTOR': 'dactolisib PI103 omipalisib',
 'mTORcat': 'AZD8055 KU0063794 torin1 torkinib',
 'rapalog': 'rapamycin temsirolimus everolimus',
 'AKT': 'MK2206 GSK690693 capivasertib ipatasertib akti_akt12 aktVIII AT13148',
 'PDPK1': 'GSK2334470 BX912', 'S6K': 'PF4708671 LY2584702 DG2',
 'SRC_ABL': 'dasatinib bosutinib imatinib', 'FGFR': 'AZD4547', 'MET_ALK': 'crizotinib PF04217903 JNJ38877605 foretinib lorlatinib ceritinib TAE684',
 'JAK': 'tofacitinib AZD1480', 'p38': 'PH797804 TAK715 SB203580 BIRB796 SB202190 SB239063', 'JNK': 'SP600125 JNKIN8 tanzisertib',
 'PKC': 'BIM1 Go6983 Go6976 enzastaurin', 'GSK3': 'tideglusib CHIR99021 LY2090314 GSK3i_BIO SB216763',
 'CDK4_6': 'ribociclib abemaciclib', 'CDK': 'AZD5438 RO3306 flavopiridol THZ531 CDK12IN3 kenpaullone cortistatinA',
 'DDR': 'VE821 AZ20 ceralasertib NU7026 NU7441 KU60019', 'CK2': 'CX4945 TBCA', 'ROCK': 'ripasudil belumosudil H1152 Y27632',
 'PAK': 'PF3758309 FRAX486', 'CAMK2': 'KN62 KN93', 'TAK1': '7Oxo', 'RTKmulti': 'nintedanib amuvatinib quizartinib pexidartinib linsitinib',
 'PKA': 'H89', 'IKK': 'SC514', 'AURK_PLK': 'MLN8054 AZD1152 VX680 BI2536 BI4834', 'RSK': 'BID1870', 'other': 'LRRK2IN1 ibrutinib',
}
A2C = {a: c for c, s in CLASS.items() for a in s.split()}
# drug binding mode (judgement from literature; flagged as such)
MODE = {'MK2206': 'allosteric', 'akti_akt12': 'allosteric', 'aktVIII': 'allosteric', 'GSK690693': 'ATP', 'capivasertib': 'ATP', 'ipatasertib': 'ATP', 'AT13148': 'ATP',
        'PD0325901': 'allosteric', 'selumetinib': 'allosteric', 'refametinib': 'allosteric', 'U0126': 'allosteric', 'trametinib': 'allosteric', 'cobimetinib': 'allosteric', 'PD184352': 'allosteric',
        'rapamycin': 'allosteric', 'temsirolimus': 'allosteric', 'everolimus': 'allosteric',
        'imatinib': 'typeII', 'sorafenib': 'typeII', 'nintedanib': 'typeII_like', 'quizartinib': 'typeII', 'lapatinib': 'typeI.5', 'vemurafenib': 'typeI.5', 'dabrafenib': 'typeI.5', 'LY3009120': 'typeII',
        'dasatinib': 'typeI', 'bosutinib': 'typeI', 'gefitinib': 'typeI', 'erlotinib': 'typeI', 'afatinib': 'covalent', 'BIRB796': 'typeII', 'tideglusib': 'non-ATP', 'TBCA': 'ATP', 'CX4945': 'ATP'}

C = load.conditions(); D = load.drug_targets().set_index('agent')
ki = C[(C.perturbation_class == 'kinase_inhibitor') & C.target_confidence.isin(['high', 'medium'])].copy()
ki = ki[~ki.qc_flag.isin(['scale_suspect', 'duplicate_vector'])]
ki = ki[ki.perturbation_subclass.isin(['inhibitor', 'inhibitor_reversed_contrast'])]
ki = ki[ki.background_agents.isna() & ~ki.is_combination.astype(bool)]
ki = ki[ki.kinase_inhibitors.isin(A2C)]
ds = ki.dose_series_id.notna() & (ki.study_pmid == '36926954')
ki = ki[~ds | ki.dose_nM.between(300, 10000)]
ki['sign'] = np.where(ki.contrast_reversed.astype(bool), -1.0, 1.0)
ki['cls'] = ki.kinase_inhibitors.map(A2C)
ki['exp'] = ki.study_pmid + '|' + ki.cell_line_norm.astype(str) + '|' + ki.kinase_inhibitors
def tset(a, col):
    if a not in D.index or pd.isna(D.loc[a, col]): return ''
    return D.loc[a, col]
E = ki.groupby('exp').agg(study=('study_pmid', 'first'), cell_line=('cell_line_norm', 'first'), agent=('kinase_inhibitors', 'first'), cls=('cls', 'first'),
                          n_cond=('conditionID', 'size'), time_min=('time_min', 'median'), dose_nM=('dose_nM', 'median'),
                          primary=('inhibitor_primary_targets', 'first')).reset_index()
E['klaeger100'] = E.agent.map(lambda a: tset(a, 'klaeger_kd_le100nM'))
E['klaeger1000'] = E.agent.map(lambda a: tset(a, 'klaeger_kd_le1000nM'))
E['hijazi'] = E.agent.map(lambda a: tset(a, 'hijazi_targets'))
E['binders'] = [';'.join(sorted(set(filter(None, ';'.join([r.primary or '', r.klaeger1000, r.hijazi]).split(';'))))) for r in E.itertuples()]
E['mode'] = E.agent.map(MODE).fillna('ATP(assumed)')
E.to_csv(OUT + 'experiments.tsv', sep='\t', index=False)
ki[['conditionID', 'exp', 'sign', 'dose_nM', 'time_min']].to_csv(OUT + 'exp_conditions.tsv', sep='\t', index=False)
print(len(ki), 'conditions ->', len(E), 'experiments;', E.study.nunique(), 'studies')
print(E.groupby('cls').agg(nexp=('exp', 'size'), nstudy=('study', 'nunique'), ncl=('cell_line', 'nunique'), nagent=('agent', 'nunique')).sort_values('nexp', ascending=False).to_string())

def agg(W):  # W: feature x condition -> feature x exp (sign-corrected mean)
    W = W[ki.conditionID] * ki.sign.values
    return W.T.groupby(ki.exp.values).mean().T
# activity
L = load.activity_long(); L = L[L.condition.isin(ki.conditionID)]
for meth, col in [('corev3', 'z_core'), ('psp_z', 'z')]:
    l = L[(L.method == meth) & (L.n_substrates >= 5)]
    W = l.pivot(index='kinase', columns='condition', values=col).reindex(columns=ki.conditionID)
    agg(W).astype('float32').to_parquet(OUT + f'act_{meth}_exp.parquet')
    if meth == 'psp_z':
        W2 = l.pivot(index='kinase', columns='condition', values='mean_log2fc').reindex(columns=ki.conditionID)
        agg(W2).astype('float32').to_parquet(OUT + 'act_pspmean_exp.parquet')
# sites
X = load.sites_matrix(ki.conditionID.tolist())
X = X[X.notna().sum(1) >= 1]
agg(X).astype('float32').to_parquet(OUT + 'site_exp.parquet')
# per-experiment global stats for normalisation-artefact control
Xe = pd.read_parquet(OUT + 'site_exp.parquet')
st = pd.DataFrame({'median': Xe.median(), 'frac_down1': (Xe < -1).mean(), 'frac_up1': (Xe > 1).mean(), 'n_sites': Xe.notna().sum(), 'mad': (Xe - Xe.median()).abs().median()})
st.to_csv(OUT + 'exp_globalstats.tsv', sep='\t')
print(st.describe())
