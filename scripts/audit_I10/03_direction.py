"""audit_I10 step 3: direction asymmetry (down- vs up-sites).
(a) Re-tabulate b_I10 T2 (their LOSO AUROCs): median down vs up, split inhibitor vs stimulus footprints.
(b) Own LOSO re-implementation (independent code, simplified): per multi-study group, hold out each study; consensus from training
    study-units (site in >=2 units, >=75% sign agreement [100% if 2 units], |mean z|>=1, Stouffer, top-200); down/up subsets; held-out
    experiments scored by sign-aligned mean z over overlapping sites (>=5); AUROC same-group vs all other labelled experiments.
    Size-matched: the larger direction subset randomly subsampled to the smaller one's size (20 draws, mean AUROC).
(c) Tumour PI3K: PI3Kaxis down/up sites (b_I10 final footprint) vs PIK3CA hotspot / PTEN mut or deep deletion; activity = -sign-aligned
    score; per-cohort AUROC + pooled AUROC on within-cohort percentiles; size-matched down subsample; matched-random-set null.
Outputs: A3_T2_retab.tsv, A3_loso_own.tsv, A3_pi3k_tumour.tsv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score
from common import *
rng = np.random.default_rng(5)
STIM = {'EGF', 'GF_all', 'GPCR_ligand', 'JAKSTAT_cytokine', 'TNF', 'DNA_damage', 'replication_stress', 'proteotoxic', 'mito_depol', 'oxidative',
        'phosphatase_inh', 'infection', 'immune_receptor', 'mitosis_cellcycle', 'serum', 'starvation', 'FGF', 'insulin', 'metformin'}
# (a)
T2 = pd.read_csv('/workspace/results/b_I10/tables/T2_loso_auroc.tsv', sep='\t')
T2['kind'] = np.where(T2.footprint.isin(STIM), 'stimulus/stress', 'inhibitor')
b = T2.dropna(subset=['auroc_down', 'auroc_up'])
ret = [dict(subset='all with both', n=len(b), med_down=b.auroc_down.median(), med_up=b.auroc_up.median(), frac_down_gt_up=(b.auroc_down > b.auroc_up).mean())]
for k, d in b.groupby('kind'):
    ret.append(dict(subset=k, n=len(d), med_down=d.auroc_down.median(), med_up=d.auroc_up.median(), frac_down_gt_up=(d.auroc_down > d.auroc_up).mean()))
pd.DataFrame(ret).to_csv(RA + 'A3_T2_retab.tsv', sep='\t', index=False); print(pd.DataFrame(ret).round(3))
# (b)
M = pd.read_csv('/workspace/cache/b_I10/fp_meta.tsv', sep='\t', dtype={'study': str}).set_index('exp')
W = pd.read_parquet('/workspace/cache/b_I10/fp_sites.parquet')[M.index]
G = pd.read_csv('/workspace/results/b_I10/tables/T0_groups.tsv', sep='\t').set_index('group')
lab = M[M.group.isin(G.index[G.eligible])]
multi = [g for g, d in lab.groupby('group') if d.study.nunique() >= 2]
Wv = np.clip(W.values.astype(np.float32), -10, 10); col = {e: i for i, e in enumerate(W.columns)}
def consensus(exps):
    m = M.loc[exps]; U = np.column_stack([np.nanmean(Wv[:, [col[e] for e in d.index]], 1) for _, d in m.groupby('study')])
    n = np.isfinite(U).sum(1)
    with np.errstate(all='ignore'):
        mu = np.nanmean(U, 1); S = np.clip(np.nansum(U, 1) / np.sqrt(n), -8, 8); agree = (np.sign(U) == np.sign(mu)[:, None]).sum(1) / n
    ok = (n >= 2) & (np.abs(mu) >= 1) & (agree >= (1 if U.shape[1] == 2 else .75))
    i = np.where(ok)[0]; i = i[np.argsort(-np.abs(S[i]))][:200]
    return i, np.sign(S[i])
def sc(ii, sg, exps):
    Z = np.clip(Wv[np.ix_(ii, [col[e] for e in exps])], -5, 5) * sg[:, None]
    n = np.isfinite(Z).sum(0); s = np.nanmean(Z, 0); s[n < 5] = np.nan; return s
rows = []
for g in multi:
    full = consensus(lab.index[lab.group == g]); gst = set(lab[lab.group == g].study)
    for st in sorted(lab.study.unique()):
        tr = lab.index[(lab.group == g) & (lab.study != st)]
        if M.loc[tr].study.nunique() < 2: continue          # need >=2 training studies for a consensus
        ii, sg = consensus(tr) if st in gst else full       # footprint never sees the held-out study
        test = lab.index[lab.study == st].tolist(); y = (lab.loc[test].group == g).values
        rows.append(dict(group=g, heldout=st, y=y, test=test, ii=ii, sg=sg))
out = []
for g in multi:
    R = [r for r in rows if r['group'] == g]
    if not R: continue
    res = {}
    for var in ['all', 'down', 'up', 'down_sm', 'up_sm']:
        Y, S_ = [], []
        reps = 20 if var.endswith('_sm') else 1
        for rep in range(reps):
            Yr, Sr = [], []
            for r in R:
                ii, sg = r['ii'], r['sg']; dn, upi = np.where(sg < 0)[0], np.where(sg > 0)[0]
                if var == 'all': sel = np.arange(len(ii))
                elif var == 'down': sel = dn
                elif var == 'up': sel = upi
                else:
                    k = min(len(dn), len(upi))
                    if k < 5: sel = np.array([], int)
                    else: sel = rng.choice(dn, k, replace=False) if var == 'down_sm' else rng.choice(upi, k, replace=False)
                if len(sel) < 5: continue
                Yr.append(r['y']); Sr.append(sc(ii[sel], sg[sel], r['test']))
            if not Yr: break
            y_, s_ = np.concatenate(Yr), np.concatenate(Sr); ok = np.isfinite(s_)
            if y_[ok].sum() >= 2 and (~y_[ok]).sum() >= 5: Y.append(roc_auc_score(y_[ok], s_[ok]))
        res[var] = np.mean(Y) if Y else np.nan
    res['n_pos'] = int(sum(r['y'].sum() for r in R)); res['n_heldout'] = len(R)
    res['med_frac_down'] = np.median([(r['sg'] < 0).mean() for r in R])
    out.append(dict(footprint=g, kind='stimulus/stress' if g in STIM else 'inhibitor', **res))
L = pd.DataFrame(out); L.to_csv(RA + 'A3_loso_own.tsv', sep='\t', index=False); print(L.round(2).to_string())
for nm, d in [('all', L), ('inhibitor', L[L.kind == 'inhibitor']), ('stimulus', L[L.kind != 'inhibitor'])]:
    b = d.dropna(subset=['down', 'up']); bs = d.dropna(subset=['down_sm', 'up_sm'])
    print(nm, 'n', len(b), 'med down/up', round(b.down.median(), 3), round(b.up.median(), 3), 'frac down>up', round((b.down > b.up).mean(), 2),
          '| size-matched n', len(bs), round(bs.down_sm.median(), 3), round(bs.up_sm.median(), 3), round((bs.down_sm > bs.up_sm).mean(), 2))
del W, Wv
# (c) tumour PI3K
P = fp_set('PI3Kaxis'); dn, up = P[P < 0], P[P > 0]
def lab_pi3k(coh):
    m = cl.table(coh, 'mutations'); pos = m.protein_change.fillna('').str.extract(r'p\.[A-Z](\d+)')[0].astype(float)
    pik = (m.gene == 'PIK3CA') & pos.isin([345, 420, 542, 545, 546, 1043, 1047, 1049, 726, 88, 38, 111])
    pten = (m.gene == 'PTEN') & m.Variant_Classification.isin(['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'In_Frame_Del', 'In_Frame_Ins'])
    pos_ = set(m.case_id[pik | pten] + '_T')
    g = cl.matrix(coh, 'cnv_gistic')
    if 'PTEN' in g.index: pos_ |= set(g.columns[g.loc['PTEN'] == -2] + '_T')
    return pos_
rows = []; pooled = {}
for coh in ['BRCA', 'UCEC', 'HNSCC', 'LSCC', 'COAD', 'GBM']:
    Z = standardise(cohort(coh)); y = Z.columns.isin(lab_pi3k(coh))
    cov = Z.notna().mean(1); sd = cohort(coh).std(1)
    bb = (pd.qcut(cov.rank(method='first'), 5, labels=False) * 5 + pd.qcut(sd.rank(method='first'), 5, labels=False))
    pool = {k: np.where(bb.values == k)[0] for k in range(25)}
    feats = {'all': P, 'down': dn, 'up': up}
    for r_ in range(20): feats[f'down_sm{r_}'] = dn[dn.index.isin(Z.index)].sample(int(up.index.isin(Z.index).sum()), random_state=r_)
    for r_ in range(100):
        u = up[up.index.isin(Z.index)]; jj = np.unique([rng.choice(pool[k]) for k in bb[u.index].values])
        feats[f'nullup{r_}'] = pd.Series(1.0, index=Z.index[jj])
        d_ = dn[dn.index.isin(Z.index)]; jj = np.unique([rng.choice(pool[k]) for k in bb[d_.index].values])
        feats[f'nulldn{r_}'] = pd.Series(-1.0, index=Z.index[jj])
    for f, s in feats.items():
        a = -score(Z, s); ok = a.notna().values
        pooled.setdefault(f, []).append(pd.DataFrame({'x': a[ok].rank(pct=True).values, 'y': y[ok]}))
        if not f.startswith(('null', 'down_sm')): rows.append(dict(cohort=coh, feature=f, auroc=roc_auc_score(y[ok], a[ok]), n_pos=int(y[ok].sum()), n_sites=int(s.index.isin(Z.index).sum())))
    print(coh, int(y.sum()), flush=True)
pa = {f: roc_auc_score(pd.concat(v).y, pd.concat(v).x) for f, v in pooled.items()}
for f in ['all', 'down', 'up']: rows.append(dict(cohort='POOLED', feature=f, auroc=pa[f]))
rows.append(dict(cohort='POOLED', feature='down size-matched to up (mean of 20)', auroc=np.mean([pa[f'down_sm{r}'] for r in range(20)])))
nu = np.array([pa[f'nullup{r}'] for r in range(100)]); nd = np.array([pa[f'nulldn{r}'] for r in range(100)])
rows.append(dict(cohort='POOLED', feature='matched random null (up-sized) median [5-95%]', auroc=np.median(nu), note=f'{np.quantile(nu,.05):.3f}-{np.quantile(nu,.95):.3f}'))
rows.append(dict(cohort='POOLED', feature='matched random null (down-sized) median [5-95%]', auroc=np.median(nd), note=f'{np.quantile(nd,.05):.3f}-{np.quantile(nd,.95):.3f}'))
rows.append(dict(cohort='POOLED', feature='p(down) vs down-sized null', auroc=(1 + (nd >= pa['down']).sum()) / 101))
rows.append(dict(cohort='POOLED', feature='p(up, two-sided) vs up-sized null', auroc=(1 + (np.abs(nu - .5) >= abs(pa['up'] - .5)).sum()) / 101))
T = pd.DataFrame(rows); T.to_csv(RA + 'A3_pi3k_tumour.tsv', sep='\t', index=False); print(T.round(3).to_string())
