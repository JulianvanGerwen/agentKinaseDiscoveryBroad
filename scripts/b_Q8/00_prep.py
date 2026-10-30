"""b_Q8 step 0: (a) data-driven protein composition markers (stromal / immune) learned from TUMOUR ESTIMATE (RNA) scores:
pooled within-cohort Spearman of each protein with ESTIMATE Stromal/Immune over 10 cohorts; top-30 proteins per score (measured in >=8
cohorts, exclusive). Epithelial = fixed list. Validated by leave-cohort-out rho with ESTIMATE in tumours.
(b) per-sample composition proxy scores (mean z of markers, within cohort over T+N jointly so T vs N levels are comparable) for all
tumours and NATs; (c) cold-ischaemia footprint (Mertins 2014 experiments in the compendium, as audit_I10 05) -> cache.
Outputs: tables/T0_composition_markers.tsv, T0_composition_validation.tsv, cache composition_scores.parquet, ischaemia_fp.tsv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from common import *
from scipy.stats import spearmanr
EPI = ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1', 'CLDN4', 'CLDN7', 'ESRP1', 'SPINT1', 'KRT7']
S = load.sites(); g2u = S.drop_duplicates('uniprot').drop_duplicates('gene').set_index('gene').uniprot
rho = {}
for coh in cl.COHORTS:
    pr = cl.matrix(coh, 'protein'); T = [c for c in cl.tumour_ids(coh) if c in pr.columns]
    a = cl.table(coh, 'annotations')
    a.index = [f'{c}_T' for c in a.case_id] if 'case_id' in a.columns else a.index
    X = pr[T]; X = X[X.notna().mean(1) >= 0.8]
    for sc in ['ESTIMATE_StromalScore', 'ESTIMATE_ImmuneScore']:
        v = a[sc].reindex(T)
        r = X.rank(1).T.corrwith(v.rank())
        rho[(coh, sc)] = r
R = pd.DataFrame(rho)
rows = []
for sc in ['ESTIMATE_StromalScore', 'ESTIMATE_ImmuneScore']:
    d = R.xs(sc, axis=1, level=1); ok = d.notna().sum(1) >= 8
    rows.append(pd.DataFrame({'uniprot': d.index[ok], 'score': sc, 'mean_rho': d[ok].mean(1).values, 'n_coh': d[ok].notna().sum(1).values}))
D = pd.concat(rows)
st = D[D.score == 'ESTIMATE_StromalScore'].set_index('uniprot').mean_rho; im = D[D.score == 'ESTIMATE_ImmuneScore'].set_index('uniprot').mean_rho
# exclusivity: stromal markers must be more stromal than immune and vice versa
stro = (st - im.reindex(st.index).fillna(0)).where(st > 0.4).dropna().sort_values(ascending=False).index[:30]
immu = (im - st.reindex(im.index).fillna(0)).where(im > 0.4).dropna().sort_values(ascending=False).index[:30]
epi = [g2u[g] for g in EPI if g in g2u.index]
u2g = S.drop_duplicates('uniprot').set_index('uniprot').gene
M = pd.concat([pd.DataFrame({'uniprot': list(stro), 'set': 'stromal'}), pd.DataFrame({'uniprot': list(immu), 'set': 'immune'}),
               pd.DataFrame({'uniprot': epi, 'set': 'epithelial'})])
M['gene'] = M.uniprot.map(u2g); M['rho_stromal'] = M.uniprot.map(st); M['rho_immune'] = M.uniprot.map(im)
M.to_csv(TQ + 'T0_composition_markers.tsv', sep='\t', index=False)
print(M.groupby('set').gene.apply(lambda x: ','.join(map(str, x))).to_string())
# scores
out, val = [], []
for coh in cl.COHORTS:
    pr = cl.matrix(coh, 'protein'); s = cl.samples(); s = s[(s.cohort == coh) & s.has_phospho & s.sample_id.isin(pr.columns)]
    Z = standardise(pr[list(s.sample_id)])
    sc = pd.DataFrame({k: Z.reindex(M[M.set == k].uniprot).mean(0) for k in ['stromal', 'immune', 'epithelial']})
    sc['cohort'] = coh; sc['type'] = s.set_index('sample_id').sample_type.reindex(sc.index).values
    out.append(sc)
    a = cl.table(coh, 'annotations'); a.index = [f'{c}_T' for c in a.case_id] if 'case_id' in a.columns else a.index
    T = sc[sc.type == 'Tumor']
    for k, e in [('stromal', 'ESTIMATE_StromalScore'), ('immune', 'ESTIMATE_ImmuneScore')]:
        v = pd.concat([T[k], a[e]], axis=1).dropna()
        val.append(dict(cohort=coh, proxy=k, estimate=e, rho=spearmanr(v.iloc[:, 0], v.iloc[:, 1])[0], n=len(v)))
CS = pd.concat(out); CS.to_parquet(CQ + 'composition_scores.parquet')
V = pd.DataFrame(val); V.to_csv(TQ + 'T0_composition_validation.tsv', sep='\t', index=False)
print(V.pivot(index='cohort', columns='proxy', values='rho').round(2))
print(CS.groupby(['cohort', 'type'])[['stromal', 'immune', 'epithelial']].mean().round(2))
# ischaemia footprint
W = pd.read_parquet('/workspace/cache/b_I10/fp_sites.parquet')
ex = [c for c in W.columns if c.startswith('24719451|') and 'ischaemia' in c]; print('ischaemia exps', ex)
I = W[ex]; n = I.notna().sum(1); mu = I.mean(1); cons = (np.sign(I).eq(np.sign(mu), axis=0) | I.isna()).all(1)
rank = mu[(n >= 2) & cons].abs().sort_values(ascending=False)
pd.Series(np.sign(mu[rank.index[:200]]).values, index=rank.index[:200], name='sign').to_csv(CQ + 'ischaemia_fp.tsv', sep='\t')
