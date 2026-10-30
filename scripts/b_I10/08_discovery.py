"""I10 step 8 (task C): tumours that mimic perturbations; reproducible genomic/clinical correlates; survival (secondary).
Score = footprint mimicry (tc_padj, null-normalised nz; + = tumour resembles the perturbed state).
Associations per cohort: OLS score ~ feature + purity (WES purity, median-imputed). Features: driver mutations (gene mutated in >=5 tumours
of the cohort and in >=3 cohorts overall), TMB(log), CIN index, ESTIMATE immune score, stage (III/IV), age.
Meta: Stouffer z (cohort-size weighted sqrt n); reproducible = BH q<0.05 and same sign with p<0.05 in >=3 cohorts.
Kinase comparison: same pipeline on core-v3 z_core (tc_padj; kinases with >=5 substrates in >=80% tumours).
Survival: Cox OS ~ score + age + late stage + purity per cohort (cohorts with >=15 events), DL random-effects meta.
Outputs: T9_assoc_footprint.tsv, T9_assoc_kinase.tsv, T10_survival.tsv"""
import sys, numpy as np, pandas as pd, warnings
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
sys.path.insert(0, '/workspace/scripts/b_I19'); from common import dl_meta
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
S = pd.read_parquet(CA + 'cptac_scores_tc_padj.parquet')
S = S[~S.footprint.str.contains(r'\|')]
F = S.pivot_table(index='sample', columns='footprint', values='nz')
KZ = cl.activity('tc_padj', 'corev3_zcore'); KN = cl.activity('tc_padj', 'corev3_n')
tum = cl.tumour_ids()
kok = ((KN[tum] >= 5).mean(1) >= 0.8); KZ = KZ.loc[kok, tum].T
print('footprints', F.shape[1], 'kinases', KZ.shape[1])
clin = cl.clinical().assign(sample=lambda d: d.case_id + '_T').set_index('sample')
ANN = pd.concat([cl.table(c, 'annotations').assign(sample=lambda d: d.case_id + '_T') for c in cl.COHORTS]).set_index('sample')
# features
feat = {}; mutc = {}
for coh in cl.COHORTS:
    mb = cl.matrix(coh, 'mut_binary'); mb.columns = mb.columns + '_T'
    ids = [s for s in cl.tumour_ids(coh) if s in mb.columns]
    m = mb[ids]; fr = m.sum(1); mutc[coh] = fr[fr >= 5].index
    feat[coh] = m
cnt = pd.Series([g for c in mutc.values() for g in c]).value_counts(); genes = cnt[cnt >= 3].index.tolist()
print('driver-like genes', len(genes))
def features(coh):
    ids = cl.tumour_ids(coh); D = pd.DataFrame(index=ids)
    m = feat[coh]
    for g in genes:
        if g in mutc[coh]: D['mut:' + g] = m.loc[g].reindex(ids).astype(float)
    c = clin.reindex(ids)
    D['TMB_log'] = np.log1p(c.TMB.astype(float)); D['CIN'] = c.CNV_index_CIN_index.astype(float)
    D['immune'] = ANN.reindex(ids).ESTIMATE_ImmuneScore.astype(float)
    D['late_stage'] = c.Stage.map({'Stage I': 0, 'Stage II': 0, 'Stage III': 1, 'Stage IV': 1}).astype(float)
    D['age'] = c.Age.astype(float)
    pur = c.WES_purity.astype(float); D['_purity'] = pur.fillna(pur.median())
    return D
def assoc(Y, tag):
    rows = []
    for coh in cl.COHORTS:
        D = features(coh); ids = D.index
        for fcol in [c for c in D.columns if not c.startswith('_')]:
            x = D[fcol]
            if x.notna().sum() < 30 or x.nunique() < 2: continue
            if fcol.startswith('mut:') and x.sum() < 5: continue
            for ycol in Y.columns:
                y = Y[ycol].reindex(ids)
                ok = x.notna() & y.notna() & D._purity.notna()
                if ok.sum() < 30: continue
                X = sm.add_constant(pd.DataFrame({'x': (x[ok] - x[ok].mean()) / (x[ok].std() if not fcol.startswith('mut:') else 1), 'pur': D._purity[ok]}))
                r = sm.OLS((y[ok] - y[ok].mean()) / y[ok].std(), X).fit()
                rows.append((ycol, fcol, coh, r.params['x'], r.tvalues['x'], r.pvalues['x'], int(ok.sum())))
    A = pd.DataFrame(rows, columns=['score', 'feature', 'cohort', 'beta', 't', 'p', 'n'])
    A['z'] = norm.isf(A.p / 2) * np.sign(A.t)
    M = A.groupby(['score', 'feature']).apply(lambda d: pd.Series({
        'n_cohorts': len(d), 'meta_z': (d.z * np.sqrt(d.n)).sum() / np.sqrt(d.n.sum()), 'mean_beta': d.beta.mean(),
        'n_pos_sig': ((d.p < 0.05) & (d.t > 0)).sum(), 'n_neg_sig': ((d.p < 0.05) & (d.t < 0)).sum(),
        'cohorts': ';'.join(f'{c}:{b:+.2f}' for c, b in zip(d.cohort, d.beta))})).reset_index()
    M['meta_p'] = 2 * norm.sf(M.meta_z.abs()); M['q'] = multipletests(M.meta_p, method='fdr_bh')[1]
    M['n_rep'] = np.where(M.meta_z > 0, M.n_pos_sig, M.n_neg_sig)
    M['reproducible'] = (M.q < 0.05) & (M.n_rep >= 3)
    M.to_csv(RT + f'T9_assoc_{tag}.tsv', sep='\t', index=False); A.to_csv(CA + f'assoc_{tag}_percohort.tsv', sep='\t', index=False)
    return M
MF = assoc(F, 'footprint'); MK = assoc(KZ, 'kinase')
# compare: per feature, best footprint vs best kinase
cmp = MF[MF.reproducible].groupby('feature').meta_z.apply(lambda z: z.abs().max()).rename('best_fp_absz').to_frame()
cmp = cmp.join(MK.groupby('feature').meta_z.apply(lambda z: z.abs().max()).rename('best_kin_absz'))
cmp = cmp.join(MK[MK.reproducible].groupby('feature').size().rename('n_kin_reproducible')).fillna({'n_kin_reproducible': 0})
cmp = cmp.join(MF[MF.reproducible].groupby('feature').size().rename('n_fp_reproducible'))
cmp.to_csv(RT + 'T9_feature_fp_vs_kinase.tsv', sep='\t')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 100)
print(MF[MF.reproducible].sort_values('meta_p').head(60).round(4).to_string())
print(cmp.sort_values('best_fp_absz', ascending=False).round(2).to_string())
print('reproducible: footprint', MF.reproducible.sum(), 'of', len(MF), '; kinase', MK.reproducible.sum(), 'of', len(MK))
# survival
rows = []
for coh in cl.COHORTS:
    ids = cl.tumour_ids(coh); c = clin.reindex(ids)
    T = c.OS_days.astype(float); E = c.OS_event.astype(float)
    if E.sum() < 15: continue
    late = c.Stage.map({'Stage I': 0, 'Stage II': 0, 'Stage III': 1, 'Stage IV': 1}).astype(float)
    pur = c.WES_purity.astype(float); age = c.Age.astype(float)
    for fp in F.columns:
        x = F[fp].reindex(ids); x = (x - x.mean()) / x.std()
        cov = pd.DataFrame({'x': x, 'age': age.fillna(age.median()), 'pur': pur.fillna(pur.median())})
        if late.notna().sum() > 0.5 * len(ids): cov['late'] = late.fillna(0)
        ok = cov.notna().all(1) & T.notna() & E.notna() & (T > 0)
        try:
            r = PHReg(T[ok].values, cov[ok].values, status=E[ok].values).fit()
            rows.append((fp, coh, r.params[0], r.bse[0], int(E[ok].sum()), int(ok.sum())))
        except Exception as e: pass
SV = pd.DataFrame(rows, columns=['footprint', 'cohort', 'logHR', 'se', 'events', 'n'])
SV.to_csv(CA + 'survival_percohort.tsv', sep='\t', index=False)
out = []
for fp, d in SV.groupby('footprint'):
    est, se, p, tau2, I2 = dl_meta(d.logHR.values, d.se.values ** 2)
    out.append(dict(footprint=fp, n_cohorts=len(d), logHR=est, se=se, p=p, I2=I2, n_nominal=int((2 * norm.sf(abs(d.logHR / d.se)) < 0.05).sum()),
                    cohorts=';'.join(f'{c}:{h:+.2f}' for c, h in zip(d.cohort, d.logHR))))
SVm = pd.DataFrame(out); SVm['q'] = multipletests(SVm.p, method='fdr_bh')[1]
SVm.sort_values('p').to_csv(RT + 'T10_survival.tsv', sep='\t', index=False)
print(SVm.sort_values('p').head(15).round(3).to_string())
