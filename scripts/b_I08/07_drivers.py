#!/usr/bin/env python3
"""I08 step 7 (Q3): what explains activity of abundance-uncoupled but variable kinases?
Targets: meanz_tc_padj kinases not classed 'abundance'/'abundance+PTM', median SNR (vs random substrate sets) >= 1.5.
Features (per cohort, tumours): other kinases' activity (meanz_tc_padj; excluding same Manning family and kinases whose
curated substrate set overlaps the target's by Jaccard>0.1), driver mutations (>=5% in cohort), CNV log2 and protein of
a curated regulator list (cyclins/CDKIs, PKA R subunits, AMPK subunits, CK2B, PP1/PP2A subunits, PTEN, DUSPs,...).
Univariate: rank partial r adjusting purity/ESTIMATE/proliferation -> DL meta -> BH FDR (per target);
reproducible = FDR<0.05, LOCO >= 80%, >=3 cohorts nominal, same sign with corev3 target activity.
Multivariate: elastic net trained on 9 cohorts, tested on the held-out cohort (Spearman), per feature block.
Outputs: tables/drivers_univariate.tsv, drivers_reproducible.tsv, drivers_enet_loco.tsv"""
import os, sys, numpy as np, pandas as pd, warnings
os.environ['OMP_NUM_THREADS'] = '8'; warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *
from scipy import stats
from sklearn.linear_model import ElasticNetCV

kin = pd.read_csv(f'{R}/model_corev3/kinases.tsv', sep='\t'); fam = dict(zip(kin.kinase, kin.manning_family)); k2acc = dict(zip(kin.kinase, kin.uniprot))
REG = ('CCNA2 CCNB1 CCNB2 CCND1 CCND2 CCND3 CCNE1 CCNE2 CCNH CCNK CCNT1 CDKN1A CDKN1B CDKN2A CDKN2B CDKN2C CDKN3 CKS1B CKS2 RB1 '
       'CDC25A CDC25B CDC25C CDK5R1 PRKAR1A PRKAR1B PRKAR2A PRKAR2B PRKAB1 PRKAB2 PRKAG1 PRKAG2 STK11 CAB39 STRADA CSNK2B '
       'PPP2CA PPP2CB PPP2R1A PPP2R2A PPP2R5A PPP2R5C PPP2R5D PPP2R5E PPP1CA PPP1CB PPP1CC PPP1R12A PPP6C PPP4C PPP3CA PPM1D '
       'PTEN PHLPP1 PHLPP2 DUSP1 DUSP4 DUSP6 PTPN11 PTPN1 PIK3CA PIK3R1 PIK3R2 KRAS NRAS HRAS NF1 MYC RPTOR RICTOR MLST8 TSC1 TSC2 '
       'DEPTOR TPX2 BORA CALM1 CALM2 CALM3 ERBB2 EGFR MET YWHAZ YWHAB SFN').split()
DRV = ('TP53 KRAS PIK3CA PTEN APC EGFR BRAF NRAS IDH1 NF1 RB1 CDKN2A STK11 KEAP1 ARID1A CTNNB1 FBXW7 PIK3R1 SMAD4 KMT2D KMT2C '
       'NFE2L2 FAT1 NOTCH1 VHL PBRM1 SETD2 BAP1 ATRX CDH1 GATA3 MAP3K1 CASP8 HRAS ATM BRCA2 PPP2R1A RNF43 TSC2 PTPN11 ERBB2 CTCF ARID1B').split()
ks = pd.read_csv(f'{R}/edges/edges_psp_only.tsv.gz', sep='\t'); ks = ks[ks.source.eq('cur')]
subs = ks.groupby('kinase').site.apply(set).to_dict()
# known upstream: curated kinase j -> site on protein of k
up = {}
for s, j in zip(ks.site, ks.kinase):
    up.setdefault(s.split('_')[0], set()).add(j)

C = pd.read_csv(f'{OUT}/tables/kinase_classification_all_sources.tsv', sep='\t')
Cm = C[C.source == 'meanz_tc_padj'].set_index('kinase')
targets = sorted(Cm[(~Cm.cls.str.startswith('abundance')) & (Cm.snr >= 1.5)].index)
print('targets', len(targets), targets, flush=True)
AM = pd.read_parquet(f'{OUT}/cache/act_meanz_tc_padj.parquet'); AC = pd.read_parquet(f'{OUT}/cache/act_corev3_tc_padj.parquet')


def jacc(a, b):
    a, b = subs.get(a, set()), subs.get(b, set()); return len(a & b) / max(1, len(a | b))


# build per-cohort feature frames (ranked, residualised) once
data = {}
for c in cl.COHORTS:
    T = [t for t in cl.tumour_ids(c) if t in AM.columns]; cases = [t[:-2] for t in T]
    Cv = covariates(c, T); Cr = np.column_stack([rank(Cv[x]) for x in Cv])
    F = {}
    A = AM[T]
    for k in A.index:
        if A.loc[k].notna().sum() >= 20: F[f'act:{k}'] = A.loc[k].values.astype(float)
    mut = cl.matrix(c, 'mut_binary').reindex(index=DRV, columns=cases)
    for g in DRV:
        if g in mut.index:
            v = mut.loc[g].values.astype(float)
            if np.nanmean(v) >= .05 and np.nansum(v) >= 5: F[f'mut:{g}'] = v
    cnv = cl.matrix(c, 'cnv_log2'); cnv = cnv[~cnv.index.duplicated()].reindex(index=REG, columns=cases)
    for g in REG:
        if g in cnv.index and cnv.loc[g].notna().sum() >= 20: F[f'cnv:{g}'] = cnv.loc[g].values.astype(float)
    prot = cl.matrix(c, 'protein'); sm = cl.site_map()[['gene', 'uniprot']].drop_duplicates('gene').set_index('gene').uniprot
    for g in REG:
        a = sm.get(g)
        if a in prot.index and prot.loc[a, T].notna().sum() >= 20: F[f'prot:{g}'] = prot.loc[a, T].values.astype(float)
    Fr = {f: resid(rank(v), Cr) for f, v in F.items()}
    tm = {}
    for src, AA in [('meanz', AM), ('corev3', AC)]:
        for k in targets:
            if k in AA.index: tm[(src, k)] = resid(rank(AA.loc[k, T].values.astype(float)), Cr)
    data[c] = dict(F=Fr, tgt=tm, n=len(T))
    print(c, len(Fr), flush=True)

rows = []
for k in targets:
    for src in ['meanz', 'corev3']:
        per = {}
        for c, d in data.items():
            y = d['tgt'].get((src, k))
            if y is None or np.isfinite(y).sum() < 20: continue
            for f, x in d['F'].items():
                if f.startswith('act:'):
                    j = f[4:]
                    if j == k or fam.get(j) == fam.get(k) or jacc(j, k) > .1: continue
                r, n = pcor(y, x)
                if np.isfinite(r): per.setdefault(f, []).append((c, r, n))
        for f, L in per.items():
            if len(L) < 3: continue
            cs, rr, nn = zip(*L); m = meta(rr, nn, 3)
            loco = sum((lambda mm: mm['p'] < .05 and np.sign(mm['r']) == np.sign(m['r']))(meta([r for cc, r in zip(cs, rr) if cc != c0], [n for cc, n in zip(cs, nn) if cc != c0], 3)) for c0 in cs)
            z = np.arctanh(np.array(rr)) * np.sqrt(np.array(nn) - 6)
            nrep = int(((np.abs(z) > 1.96) & (np.sign(rr) == np.sign(m['r']))).sum())
            rows.append(dict(target=k, src=src, feature=f, block=f.split(':')[0], r=m['r'], lo=m['lo'], hi=m['hi'], p=m['p'], I2=m['I2'],
                             k=m['k'], loco=loco / m['k'], n_rep=nrep,
                             known_upstream=(f.startswith('act:') and f[4:] in up.get(k2acc.get(k, ''), set()))))
U = pd.DataFrame(rows)
U['fdr'] = np.nan
for (t, s), idx in U.groupby(['target', 'src']).groups.items(): U.loc[idx, 'fdr'] = bh(U.loc[idx, 'p'])
U.to_csv(f'{OUT}/tables/drivers_univariate.tsv', sep='\t', index=False)
W = U.pivot_table(index=['target', 'feature', 'block'], columns='src', values=['r', 'lo', 'hi', 'fdr', 'loco', 'n_rep', 'k'], aggfunc='first')
W.columns = [f'{a}_{b}' for a, b in W.columns]; W = W.reset_index()
W['known_upstream'] = W.apply(lambda r: r.feature.startswith('act:') and r.feature[4:] in up.get(k2acc.get(r.target, ''), set()), axis=1)
W['reproducible'] = (W.fdr_meanz < .05) & (W.loco_meanz >= .8) & (W.n_rep_meanz >= 3) & (np.sign(W.r_meanz) == np.sign(W.r_corev3)) & (W.fdr_corev3 < .05)
W.sort_values(['reproducible', 'fdr_meanz'], ascending=[False, True]).to_csv(f'{OUT}/tables/drivers_reproducible.tsv', sep='\t', index=False)

# elastic net LOCO per block
blocks = {'kinase_act': ['act'], 'mutation': ['mut'], 'cnv_reg': ['cnv'], 'prot_reg': ['prot'], 'all': ['act', 'mut', 'cnv', 'prot']}
er = []
for k in targets:
    feats = sorted({f for d in data.values() for f in d['F']
                    if not (f.startswith('act:') and (f[4:] == k or fam.get(f[4:]) == fam.get(k) or jacc(f[4:], k) > .1))})
    Xc, yc = {}, {}
    for c, d in data.items():
        y = d['tgt'].get(('meanz', k))
        if y is None: continue
        X = np.column_stack([d['F'].get(f, np.full(d['n'], np.nan)) for f in feats])
        X = (X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9); X = np.nan_to_num(X)
        m = np.isfinite(y)
        if m.sum() < 20: continue
        Xc[c] = X[m]; yc[c] = (y[m] - y[m].mean()) / (y[m].std() + 1e-9)
    for b, pre in blocks.items():
        cols = [i for i, f in enumerate(feats) if f.split(':')[0] in pre]
        if not cols: continue
        for c in Xc:
            tr = [x for x in Xc if x != c]
            Xtr = np.vstack([Xc[x][:, cols] for x in tr]); ytr = np.concatenate([yc[x] for x in tr])
            mdl = ElasticNetCV(l1_ratio=.5, alphas=15, cv=3, max_iter=3000, n_jobs=1).fit(Xtr, ytr)
            pr = mdl.predict(Xc[c][:, cols])
            rho = stats.spearmanr(pr, yc[c])[0] if np.std(pr) > 0 else 0.0
            top = [feats[cols[i]] for i in np.argsort(-np.abs(mdl.coef_))[:5] if mdl.coef_[i] != 0] if b == 'all' else []
            er.append(dict(target=k, block=b, heldout=c, rho=rho, n_nonzero=int((mdl.coef_ != 0).sum()), top='; '.join(top)))
    print('enet', k, flush=True)
E = pd.DataFrame(er); E.to_csv(f'{OUT}/tables/drivers_enet_loco.tsv', sep='\t', index=False)
print(E.groupby(['target', 'block']).rho.median().unstack().round(2).to_string())
pd.set_option('display.width', 250)
print(W[W.reproducible].sort_values('fdr_meanz').groupby('target').head(4)[['target', 'feature', 'r_meanz', 'lo_meanz', 'hi_meanz', 'fdr_meanz', 'n_rep_meanz', 'k_meanz', 'r_corev3', 'known_upstream']].round(3).to_string())
