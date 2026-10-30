"""I10 step 4: validate CPTAC footprint scores against genomic/immune ground truth, compare with kinase activity (core-v3).
Pathway activity = +score for stimulus footprints (EGF, GF_all, insulin, TNF, JAKSTAT), -score for inhibitor footprints.
Per cohort AUROC (binary) or Spearman (continuous); pooled AUROC on within-cohort percentiles; null = 200 matched random site sets
(same weights; pooled null computed on the same random-set index across cohorts). Added value: logistic LRT with cohort fixed effects,
y ~ cohort + kinase-activity (+ footprint) and vice versa.  Output: T4_validation.tsv, T4_added_value.tsv"""
import sys, numpy as np, pandas as pd, warnings, re
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr, chi2
import statsmodels.api as sm
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
names = pd.read_csv(CA + 'cptac_score_names.tsv', sep='\t', header=None)[0].tolist()
STIM = ('EGF', 'GF_all', 'insulin', 'TNF', 'JAKSTAT_cytokine', 'FGF', 'serum')
sgn = lambda f: 1 if f.split('|')[0] in STIM else -1

def hot(coh):
    m = cl.table(coh, 'mutations'); m = m[m.Variant_Classification.isin(['Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site'])]
    pc = m.protein_change.fillna('')
    pos = pc.str.extract(r'p\.[A-Z](\d+)')[0].astype(float)
    ras = ((m.gene.isin(['KRAS', 'NRAS', 'HRAS'])) & pos.isin([12, 13, 61, 117, 146])) | ((m.gene == 'BRAF') & pos.isin([600, 601, 469, 464, 466, 597]))
    pik = (m.gene == 'PIK3CA') & pos.isin([345, 420, 542, 545, 546, 1043, 1047, 1049, 726, 88, 38, 111])
    pten = (m.gene == 'PTEN')
    egfrm = (m.gene == 'EGFR') & (m.Variant_Classification.isin(['Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins']))
    g = cl.matrix(coh, 'cnv_gistic')
    lab = pd.DataFrame(index=g.columns)
    lab['RAS_RAF'] = lab.index.isin(m.case_id[ras])
    ptendel = g.loc['PTEN'] == -2 if 'PTEN' in g.index else False
    lab['PI3K_PTEN'] = lab.index.isin(m.case_id[pik | pten]) | ptendel
    lab['EGFR_amp'] = (g.loc['EGFR'] == 2).values if 'EGFR' in g.index else False
    lab['EGFR_amp_or_mut'] = lab.EGFR_amp | lab.index.isin(m.case_id[egfrm])
    lab.index = lab.index + '_T'; return lab

S = pd.concat([pd.read_parquet(CA + 'cptac_scores_tc_padj.parquet').assign(variant='tc_padj'),
               pd.read_parquet(CA + 'cptac_scores_tc.parquet').assign(variant='tc')])
S['act'] = S.nz * S.footprint.map(sgn)
KIN = {'MAPK1/3': ['MAPK1', 'MAPK3'], 'EGFR': ['EGFR'], 'AKT': ['AKT1', 'AKT2'], 'MTOR': ['MTOR'], 'JAK': ['JAK1', 'JAK2', 'TYK2'],
       'IKBKB/TBK1': ['IKBKB', 'TBK1'], 'MAP2K1': ['MAP2K1'], 'RPS6KB1': ['RPS6KB1']}
KA = {v: cl.activity(v, 'corev3_z') for v in ['tc_padj', 'tc']}
kact = lambda v, k: KA[v].loc[KIN[k]].mean(0)
LAB = pd.concat([hot(c).assign(cohort=c) for c in cl.COHORTS])
ANN = pd.concat([cl.table(c, 'annotations').assign(sample=lambda d: d.case_id + '_T') for c in cl.COHORTS]).set_index('sample')

TESTS = [  # (name, label/annotation, cohorts, footprints, kinase baselines, kind)
 ('RAS/RAF hotspot', 'RAS_RAF', ['LUAD', 'COAD', 'UCEC', 'PDAC', 'LSCC', 'HNSCC'], ['KI_MEK', 'KI_MEK|down', 'KI_MEK|up', 'MAPKaxis', 'MAPKaxis|down', 'KI_ERK', 'KI_RAF', 'KI_EGFR', 'EGF', 'GF_all'], ['MAPK1/3', 'MAP2K1'], 'bin'),
 ('EGFR amp (GBM)', 'EGFR_amp', ['GBM'], ['EGF', 'GF_all', 'KI_EGFR', 'KI_EGFR|down', 'KI_MEK', 'MAPKaxis'], ['EGFR', 'MAPK1/3'], 'bin'),
 ('EGFR amp/mut (LUAD)', 'EGFR_amp_or_mut', ['LUAD'], ['EGF', 'GF_all', 'KI_EGFR', 'KI_EGFR|down', 'KI_MEK'], ['EGFR', 'MAPK1/3'], 'bin'),
 ('PIK3CA hotspot/PTEN', 'PI3K_PTEN', ['BRCA', 'UCEC', 'HNSCC', 'LSCC', 'COAD', 'GBM'], ['PI3Kaxis', 'PI3Kaxis|down', 'PI3Kaxis|up', 'KI_PI3K', 'KI_AKT', 'KI_AKT|down', 'KI_mTORcat', 'KI_PI3K_mTOR', 'insulin', 'GF_all'], ['AKT', 'MTOR', 'RPS6KB1'], 'bin'),
 ('Immune (ESTIMATE)', 'ESTIMATE_ImmuneScore', None, ['JAKSTAT_cytokine', 'TNF', 'KI_JAK', 'immune_receptor', 'infection'], ['JAK', 'IKBKB/TBK1'], 'cont'),
 ('IFNg hallmark', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', None, ['JAKSTAT_cytokine', 'TNF', 'KI_JAK', 'infection'], ['JAK', 'IKBKB/TBK1'], 'cont'),
 ('TNF/NFkB hallmark', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', None, ['TNF', 'JAKSTAT_cytokine', 'GF_all'], ['IKBKB/TBK1', 'JAK'], 'cont'),
 ('PROGENy MAPK (RNA)', 'PROGENy_MAPK', None, ['KI_MEK', 'MAPKaxis', 'KI_MEK|down', 'EGF', 'GF_all'], ['MAPK1/3', 'MAP2K1'], 'cont'),
 ('PROGENy EGFR (RNA)', 'PROGENy_EGFR', None, ['EGF', 'KI_EGFR', 'KI_MEK', 'GF_all'], ['EGFR', 'MAPK1/3'], 'cont'),
 ('PROGENy PI3K (RNA)', 'PROGENy_PI3K', None, ['PI3Kaxis', 'KI_AKT', 'insulin'], ['AKT', 'MTOR'], 'cont'),
]
rows, add = [], []
for tname, lab, cohs, fps, kins, kind in TESTS:
    cohs = cohs or cl.COHORTS
    for v in ['tc_padj', 'tc']:
        for feat in fps + ['kin:' + k for k in kins] + ['ann:PTM_SEA_PERT-PSP_EGF']:
            per = []; pooled = []
            for coh in cohs:
                ids = [s for s in cl.tumour_ids(coh)]
                if feat.startswith('kin:'):
                    x = kact(v, feat[4:]).reindex(ids)
                elif feat.startswith('ann:'):
                    x = ANN.reindex(ids)[feat[4:]].astype(float)
                else:
                    d = S[(S.variant == v) & (S.cohort == coh) & (S.footprint == feat)]
                    if d.empty: continue
                    x = d.set_index('sample').act.reindex(ids)
                y = (LAB.reindex(ids)[lab] if kind == 'bin' else ANN.reindex(ids)[lab]).astype(float)
                ok = x.notna() & y.notna(); x, y = x[ok], y[ok]
                if kind == 'bin':
                    if y.sum() < 4 or (1 - y).sum() < 4: continue
                    e = roc_auc_score(y, x)
                else:
                    if len(y) < 20: continue
                    e = spearmanr(x, y)[0]
                per.append((coh, e, int(y.sum()) if kind == 'bin' else len(y)))
                pooled.append(pd.DataFrame({'x': x.rank(pct=True), 'y': y, 'coh': coh}))
            if not per: continue
            P = pd.concat(pooled)
            pe = roc_auc_score(P.y, P.x) if kind == 'bin' else np.mean([p[1] for p in per])
            # null for footprints: same random set index pooled across cohorts
            pnull = np.nan
            if not feat.startswith(('kin:', 'ann:')):
                nulls = []
                j = names.index(feat)
                for coh in [p[0] for p in per]:
                    arr = np.load(CA + f'cptac_null_{v}_{coh}.npy')[j]  # NR x samples
                    ids = cl.tumour_ids(coh)
                    tum = pd.read_parquet(CA + f'cptac_rank_{v}_{coh}.parquet', columns=None).columns
                    nd = pd.DataFrame(arr.T * sgn(feat), index=tum)
                    nulls.append(nd)
                N = pd.concat(nulls)
                yy = (LAB.reindex(N.index)[lab] if kind == 'bin' else ANN.reindex(N.index)[lab]).astype(float)
                coh_of = pd.concat([pd.Series(p[0], index=nl.index) for p, nl in zip(per, nulls)])
                ok = yy.notna()
                ne = []
                for k in range(N.shape[1]):
                    xk = N[k][ok].groupby(coh_of[ok]).rank(pct=True)
                    if kind == 'bin': ne.append(roc_auc_score(yy[ok], xk.fillna(0.5)))
                    else: ne.append(np.nanmean([spearmanr(N[k][ok & (coh_of == c)], yy[ok & (coh_of == c)], nan_policy='omit')[0] for c in coh_of.unique()]))
                ne = np.array(ne); pnull = (1 + (ne >= pe).sum()) / (1 + len(ne))
            rows.append(dict(test=tname, variant=v, feature=feat, pooled=pe, p_null=pnull, n_cohorts=len(per),
                             per_cohort=';'.join(f'{c}:{e:.2f}' for c, e, n in per), n_pos=sum(p[2] for p in per)))
        # added value (binary tests): logistic LRT with cohort fixed effects
        if kind == 'bin':
            for fp in fps[:4]:
                for k in kins[:1]:
                    D = []
                    for coh in cohs:
                        ids = cl.tumour_ids(coh)
                        d = S[(S.variant == v) & (S.cohort == coh) & (S.footprint == fp)].set_index('sample').act.reindex(ids)
                        D.append(pd.DataFrame({'fp': d, 'kin': kact(v, k).reindex(ids), 'y': LAB.reindex(ids)[lab].astype(float), 'coh': coh}))
                    D = pd.concat(D).dropna()
                    D = D[D.groupby('coh').y.transform(lambda s: (s.sum() >= 4) & ((1 - s).sum() >= 4))]
                    if D.empty: continue
                    for c in ['fp', 'kin']: D[c] = D.groupby('coh')[c].transform(lambda s: (s - s.mean()) / s.std())
                    Xc = pd.get_dummies(D.coh, drop_first=True, dtype=float)
                    def ll(cols):
                        Xm = sm.add_constant(pd.concat([Xc, D[cols]], axis=1) if cols else Xc, has_constant='add')
                        return sm.Logit(D.y, Xm).fit(disp=0).llf
                    l0, lk, lf, lb = ll([]), ll(['kin']), ll(['fp']), ll(['kin', 'fp'])
                    add.append(dict(test=tname, variant=v, footprint=fp, kinase=k, n=len(D), n_pos=int(D.y.sum()),
                                    p_fp_beyond_kin=chi2.sf(2 * (lb - lk), 1), p_kin_beyond_fp=chi2.sf(2 * (lb - lf), 1),
                                    dev_kin=2 * (lk - l0), dev_fp=2 * (lf - l0), dev_both=2 * (lb - l0)))
    print(tname, flush=True)
R = pd.DataFrame(rows); R.to_csv(RT + 'T4_validation.tsv', sep='\t', index=False)
A = pd.DataFrame(add); A.to_csv(RT + 'T4_added_value.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80)
print(R.round(3).to_string()); print(A.round(4).to_string())
