"""audit_I10 step 2: is the MEK/EGF '-0.88' tumour axis built in, and what is it biologically?
Per cohort (tc_padj, also tc): own footprint scores (signed mean of standardised sites).
 (i) r(MEK, EGF) full; without the 35 shared sites; without shared proteins; MEK-down-only vs EGF-up-only (all sites same sign within set).
 (ii) reference: r between two disjoint random halves of PSP MAPK1/3 substrates (all +), and between two random matched site sets.
 (iii) MAPK output = -MEK score: corr with core-v3 MAPK1/3 (tc_padj), PROGENy MAPK, RAS/RAF hotspot AUROC, sample PC1 of all sites,
       purity, ESTIMATE stromal/immune, Hallmark E2F, ischaemia proxy (mean z of immediate-early/stress phosphosites, tc).
 (iv) b_I10 cached footprint-score matrix (nz): PCA PC1 variance & loadings; PC1 vs MAPK1/3; PC1 % for 49 random-site score sets.
Outputs: A2_axis_percohort.tsv, A2_correlates.tsv, A2_pca.tsv"""
import numpy as np, pandas as pd, warnings, re
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from common import *
rng = np.random.default_rng(3)
SA = site_annot()
MEK, EGF = fp_set('KI_MEK'), fp_set('EGF')
shared = MEK.index.intersection(EGF.index)
protM = set(s.split('_')[0] for s in MEK.index); protE = set(s.split('_')[0] for s in EGF.index)
PSP = psp_sets(5); MAPK = pd.concat([PSP.get('MAPK1', pd.Series(dtype=float)), PSP.get('MAPK3', pd.Series(dtype=float))])
MAPK = MAPK[~MAPK.index.duplicated()]
KA = cl.activity('tc_padj', 'corev3_z'); KAtc = cl.activity('tc', 'corev3_z')
IEG = ['JUN', 'JUNB', 'JUND', 'FOS', 'FOSB', 'FOSL1', 'FOSL2', 'EGR1', 'ATF2', 'ATF3', 'NR4A1', 'IER2', 'IER3', 'ZFP36', 'HSPB1', 'MAPKAPK2', 'CREB1', 'SRF', 'DUSP1', 'KLF6']
iegu = set(SA[SA.gene.isin(IEG)].uniprot)
def hot(coh):
    m = cl.table(coh, 'mutations'); pos = m.protein_change.fillna('').str.extract(r'p\.[A-Z](\d+)')[0].astype(float)
    ras = (m.gene.isin(['KRAS', 'NRAS', 'HRAS']) & pos.isin([12, 13, 61, 117, 146])) | ((m.gene == 'BRAF') & pos.isin([600, 601, 469, 464, 466, 597]))
    return set(m.case_id[ras] + '_T')
rows, cor = [], []
for coh in COH:
    Zp = standardise(cohort(coh)); Zt = standardise(cohort(coh, 'input_tc'))
    sM, sE = score(Zp, MEK), score(Zp, EGF)
    r = lambda a, b: pd.concat([a, b], axis=1).dropna().corr().iloc[0, 1]
    d = dict(cohort=coh, r_full=r(sM, sE),
             r_noshared=r(score(Zp, MEK.drop(shared)), score(Zp, EGF.drop(shared))),
             r_noshared_prot=r(score(Zp, MEK[[s.split('_')[0] not in protE for s in MEK.index]]), score(Zp, EGF[[s.split('_')[0] not in protM for s in EGF.index]])),
             r_MEKdown_EGFup_noshared=r(score(Zp, MEK[MEK < 0].drop(shared, errors='ignore')), score(Zp, EGF[EGF > 0].drop(shared, errors='ignore'))),
             r_full_tc=r(score(Zt, MEK), score(Zt, EGF)))
    # PSP MAPK1/3 split halves
    m = MAPK[MAPK.index.isin(Zp.index)]; hs = []
    for _ in range(50):
        p = rng.permutation(m.index); a, b = p[:len(p) // 2], p[len(p) // 2:]
        hs.append(r(score(Zp, m[a], 5), score(Zp, m[b], 5)))
    d['r_PSP_MAPK13_halves'] = np.nanmedian(hs); d['n_PSP_MAPK13'] = len(m)
    # random matched pairs (coverage x sd bins), same signs as MEK and EGF
    cov = Zp.notna().mean(1); x = cohort(coh); sd = x.std(1)
    b = (pd.qcut(cov.rank(method='first'), 5, labels=False) * 5 + pd.qcut(sd.rank(method='first'), 5, labels=False))
    pool = {k: np.where(b.values == k)[0] for k in range(25)}
    rr = []
    for _ in range(50):
        def rs(S):
            S = S[S.index.isin(Zp.index)]
            jj = [rng.choice(pool[k]) for k in b[S.index].values]
            return pd.Series(S.values, index=Zp.index[jj])
        a_, b_ = rs(MEK), rs(EGF); rr.append(r(score(Zp, a_[~a_.index.duplicated()]), score(Zp, b_[~b_.index.duplicated()])))
    d['r_random_matched_median'] = np.nanmedian(rr); d['r_random_matched_q05_q95'] = f'{np.nanquantile(rr, .05):.2f}..{np.nanquantile(rr, .95):.2f}'
    rows.append(d)
    # correlates of MAPK output (= -MEK score)
    out = -sM
    A = Zp[Zp.notna().mean(1) >= 0.5].fillna(0); U, s, Vt = np.linalg.svd(A.values - A.values.mean(1, keepdims=True), full_matrices=False)
    pc = pd.Series(Vt[0], index=Zp.columns); pc2 = pd.Series(Vt[1], index=Zp.columns)
    ann = cl.table(coh, 'annotations').assign(s=lambda q: q.case_id + '_T').set_index('s')
    cli = cl.table(coh, 'clinical').assign(s=lambda q: q.case_id + '_T').set_index('s')
    ieg = Zt[[s.split('_')[0] in iegu for s in Zt.index]].mean(0)
    C = {'MAPK1/3 corev3 (tc_padj)': KA.loc[['MAPK1', 'MAPK3']].mean(0), 'MAPK1/3 corev3 (tc)': KAtc.loc[['MAPK1', 'MAPK3']].mean(0),
         'EGF footprint (own)': sE, 'PROGENy_MAPK': ann.PROGENy_MAPK, 'PROGENy_EGFR': ann.PROGENy_EGFR, 'all-site PC1 (|r|)': pc, 'all-site PC2 (|r|)': pc2,
         'WES_purity': cli.WES_purity, 'ESTIMATE_Stromal': ann.ESTIMATE_StromalScore, 'ESTIMATE_Immune': ann.ESTIMATE_ImmuneScore,
         'HALLMARK_E2F': ann.HALLMARK_E2F_TARGETS, 'ischaemia proxy (IEG/stress sites, tc)': ieg, 'n_ieg_sites': None}
    for k, v in C.items():
        if v is None: continue
        v = pd.to_numeric(v.reindex(out.index), errors='coerce'); ok = v.notna() & out.notna()
        rho = spearmanr(out[ok], v[ok])[0] if ok.sum() > 20 else np.nan
        if '(|r|)' in k: rho = abs(rho)
        cor.append(dict(cohort=coh, var=k, rho=rho, n=int(ok.sum())))
    hs_ = hot(coh); y = out.index.isin(hs_)
    if 4 <= y.sum() <= len(y) - 4:
        ok = out.notna().values
        cor.append(dict(cohort=coh, var='RAS/RAF hotspot AUROC (MAPK output)', rho=roc_auc_score(y[ok], out[ok]), n=int(y.sum())))
        kk = KA.loc[['MAPK1', 'MAPK3']].mean(0).reindex(out.index)
        cor.append(dict(cohort=coh, var='RAS/RAF hotspot AUROC (MAPK1/3 corev3 tc_padj)', rho=roc_auc_score(y[kk.notna().values], kk.dropna()), n=int(y.sum())))
    print(coh, {k: round(v, 2) if isinstance(v, float) else v for k, v in d.items()}, flush=True)
R = pd.DataFrame(rows); R.to_csv(RA + 'A2_axis_percohort.tsv', sep='\t', index=False)
Cr = pd.DataFrame(cor); Cr.to_csv(RA + 'A2_correlates.tsv', sep='\t', index=False)
print(R.round(2).to_string()); print(Cr.groupby('var').rho.agg(['median', 'min', 'max', 'size']).round(2).to_string())
# (iv) PCA on b_I10 cached scores
S = pd.read_parquet('/workspace/cache/b_I10/cptac_scores_tc_padj.parquet'); S = S[~S.footprint.str.contains(r'\|')]
F = S.pivot_table(index='sample', columns='footprint', values='nz')
Fz = F.groupby(S.groupby('sample').cohort.first().reindex(F.index)).transform(lambda c: (c - c.mean()) / c.std()).fillna(0)
U, s, Vt = np.linalg.svd(Fz.values, full_matrices=False); ve = s ** 2 / (s ** 2).sum()
pc1 = pd.Series(U[:, 0] * s[0], index=Fz.index); load_ = pd.Series(Vt[0], index=Fz.columns)
if load_['KI_MEK'] > 0: pc1, load_ = -pc1, -load_
mk = KA.loc[['MAPK1', 'MAPK3']].mean(0).reindex(pc1.index)
out = [dict(what='PC1 var explained (49 footprint nz scores)', value=ve[0]), dict(what='PC2 var explained', value=ve[1]),
       dict(what='rho(PC1, MAPK1/3 corev3 tc_padj), pooled within-cohort z', value=spearmanr(pc1, mk, nan_policy='omit')[0])]
# Effective redundancy: footprint-site Jaccard among top PC1 loaders
top = load_.abs().sort_values(ascending=False).index[:12]
out.append(dict(what='top-12 |PC1| loaders', value=';'.join(f'{t}:{load_[t]:+.2f}' for t in top)))
# PCA after dropping footprints sharing >10% sites with KI_MEK or EGF
sets = {g: set(FP[FP.group == g].site) for g in FP.group.unique()}
red = [g for g in Fz.columns if g in sets and (len(sets[g] & sets['KI_MEK']) > 20 or len(sets[g] & sets['EGF']) > 20)]
keep = [g for g in Fz.columns if g not in red] + ['KI_MEK']
U2, s2, V2 = np.linalg.svd(Fz[keep].values, full_matrices=False)
out.append(dict(what=f'PC1 var after collapsing MAPK-overlapping footprints to KI_MEK only (dropped {len(red)-1}: {";".join(sorted(set(red)-{"KI_MEK"}))})', value=(s2 ** 2 / (s2 ** 2).sum())[0]))
out.append(dict(what='... its top loaders', value=';'.join(f'{t}:{v:+.2f}' for t, v in pd.Series(V2[0], index=keep).sort_values(key=abs, ascending=False)[:8].items())))
P = pd.DataFrame(out); P.to_csv(RA + 'A2_pca.tsv', sep='\t', index=False); pd.set_option('display.max_colwidth', 300); print(P.to_string())
