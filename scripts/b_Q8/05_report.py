"""b_Q8 step 5: summary tables + figures for Q1 (coherence), MEK/EGF axis, Q2 (edges) and Q3 (confounds)."""
import numpy as np, pandas as pd, glob, warnings
warnings.filterwarnings('ignore')
from common import *
from scipy.stats import norm, chi2, spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = np.empty_like(p); m = np.isfinite(p).sum()
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, len(p) + 1))[::-1])[::-1]; return np.minimum(q, 1)
# ---------------- Q1 coherence ----------------
C = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(TQ + 'T1_coherence_*.tsv')]).drop_duplicates(['set', 'cohort', 'cond'])
C.to_csv(TQ + 'T1_coherence.tsv', sep='\t', index=False)
rows = []
for (s, cond), d in C.groupby(['set', 'cond']):
    m = d[d.cohort.isin(MAINC)]
    if len(m) < 4: continue
    r = dict(set=s, cond=cond, k=len(m), coh_N_med=m.coh_N.median(), coh_T_med=m.coh_T.median(), z_N_med=m.z_N.median(), z_T_med=m.z_T.median(),
             nN_sig=(m.p_N < 0.05).sum(), nT_sig=(m.p_T < 0.05).sum(), delta_med=m.delta.median(), z_delta_null_med=m.z_delta_null.median(),
             all8_nN_sig=(d.p_N < 0.05).sum(), all8_nT_sig=(d.p_T < 0.05).sum(), k8=len(d))
    if 'z_delta_perm' in m and m.z_delta_perm.notna().all():
        z = m.z_delta_perm.values; r.update(meta_z_perm=z.sum() / np.sqrt(len(z)), n_up=((m.p_delta_perm < 0.05) & (z > 0)).sum(),
                                            n_down=((m.p_delta_perm < 0.05) & (z < 0)).sum(), Q_p=chi2.sf(((z - z.mean()) ** 2).sum(), len(z) - 1))
        r['meta_p'] = 2 * norm.sf(abs(r['meta_z_perm']))
    rows.append(r)
S = pd.DataFrame(rows)
for cond in S.cond.unique():
    m = (S.cond == cond) & S.meta_p.notna(); S.loc[m, 'meta_q'] = bh(S.loc[m, 'meta_p'])
def klass(r):
    nat = r.nN_sig >= 4; tum = r.nT_sig >= 4
    if pd.notna(r.meta_q) and r.meta_q < 0.05 and r.meta_z_perm > 0 and r.n_up >= 3: g = 'GAINED in tumour'
    elif pd.notna(r.meta_q) and r.meta_q < 0.05 and r.meta_z_perm < 0 and r.n_down >= 3: g = 'LOST/attenuated in tumour'
    else: g = 'preserved' if (nat and tum) else ('tumour-only (n.s. diff)' if tum else ('NAT-only (n.s. diff)' if nat else 'incoherent'))
    return g
S['class'] = S.apply(klass, axis=1)
S.to_csv(TQ + 'T2_coherence_summary.tsv', sep='\t', index=False)
P = S[S.cond == 'padj_pc3'].set_index('set')
X = S[S.cond == 'padj_comp_pc3'].set_index('set')
P['class_comp_pc3'] = X['class']; P['meta_z_comp_pc3'] = X.meta_z_perm
P['meta_z_padj'] = S[S.cond == 'padj'].set_index('set').meta_z_perm
P['delta_med_pc5'] = S[S.cond == 'padj_pc5'].set_index('set').delta_med
P['delta_med_raw'] = S[S.cond == 'raw'].set_index('set').delta_med
P = P.sort_values('meta_z_perm')
cols = ['nN_sig', 'nT_sig', 'z_N_med', 'z_T_med', 'coh_N_med', 'coh_T_med', 'delta_med', 'meta_z_perm', 'n_up', 'n_down', 'Q_p', 'meta_q', 'class', 'meta_z_padj', 'meta_z_comp_pc3', 'class_comp_pc3', 'delta_med_pc5', 'delta_med_raw', 'all8_nN_sig', 'all8_nT_sig']
print(P[cols].round(3).to_string())
P[cols].to_csv(TQ + 'T2_coherence_classes_primary.tsv', sep='\t')
print(P['class'].value_counts(), P.class_comp_pc3.value_counts())
# tumour subsample sensitivity: paired tumours vs random tumour subsamples
cc = C[C.cond == 'padj_pc3']; print('rho(coh_T paired, coh_T random subsample)', spearmanr(cc.coh_T, cc.coh_Tsub_mean, nan_policy='omit')[0],
                                    'median diff', (cc.coh_T - cc.coh_Tsub_mean).median())
# ---------------- MEK/EGF axis ----------------
F = pd.concat([pd.read_parquet(f) for f in glob.glob(CQ + 'fp_scores_*.parquet')])
ax = []
for (coh, grp, cond), d in F.groupby(['cohort', 'grp', 'cond']):
    fp = d[[c for c in d.columns if c.startswith('FP:')]].dropna(axis=1, thresh=int(0.8 * len(d)))
    fp = fp.fillna(fp.mean()); zz = (fp - fp.mean()) / fp.std()
    U, s, Vt = np.linalg.svd(zz.values, full_matrices=False); ve = s ** 2 / (s ** 2).sum()
    load_ = pd.Series(Vt[0], index=fp.columns); load_ *= np.sign(load_.get('FP:EGF', 1))
    ax.append(dict(cohort=coh, grp=grp, cond=cond, n=len(d), r_MEK_EGF_fullsets=d['FP:KI_MEK'].corr(d['FP:EGF']), PC1_var=ve[0], PC2_var=ve[1],
                   PC1_top=';'.join(f'{k[3:]}:{v:+.2f}' for k, v in load_.abs().sort_values(ascending=False).head(6).items()),
                   PC1_MEK_load=load_.get('FP:KI_MEK'), PC1_EGF_load=load_.get('FP:EGF'), loadings=load_.to_dict()))
A = pd.DataFrame(ax)
sim = []
for (coh, cond), d in A.groupby(['cohort', 'cond']):
    if len(d) < 2: continue
    lt = pd.Series(d[d.grp == 'T'].loadings.iloc[0]); ln = pd.Series(d[d.grp == 'N'].loadings.iloc[0]); k_ = lt.index.intersection(ln.index)
    sim.append(dict(cohort=coh, cond=cond, cos_PC1_T_N=abs(np.dot(lt[k_], ln[k_]) / np.linalg.norm(lt[k_]) / np.linalg.norm(ln[k_]))))
A = A.merge(pd.DataFrame(sim), on=['cohort', 'cond']).drop(columns='loadings')
# exclusive-site MEK-EGF from edge stats
E0 = pd.read_parquet(CQ + 'edge_stats_padj_pc3.parquet'); E1 = pd.read_parquet(CQ + 'edge_stats_padj.parquet')
for nm, EE in [('padj_pc3', E0), ('padj', E1)]:
    me = EE[EE.rel == 'FP:EGF|FP:KI_MEK'].set_index('cohort')
    A.loc[(A.cond == nm) & (A.grp == 'T'), 'r_MEK_EGF_exclusive'] = A[(A.cond == nm) & (A.grp == 'T')].cohort.map(me.r_T)
    A.loc[(A.cond == nm) & (A.grp == 'N'), 'r_MEK_EGF_exclusive'] = A[(A.cond == nm) & (A.grp == 'N')].cohort.map(me.r_N)
A.to_csv(TQ + 'T4_axis.tsv', sep='\t', index=False)
print(A[A.cond.isin(['padj', 'padj_pc3'])].sort_values(['cond', 'cohort', 'grp'])[['cohort', 'grp', 'cond', 'n', 'r_MEK_EGF_fullsets', 'r_MEK_EGF_exclusive', 'PC1_var', 'PC1_MEK_load', 'cos_PC1_T_N', 'PC1_top']].round(2).to_string())
# ---------------- Q3: composition / ischaemia in NAT vs tumour ----------------
CS = pd.read_parquet(CQ + 'composition_scores.parquet')
q3 = []
for coh in NATC:
    d = F[(F.cohort == coh) & (F.cond == 'padj')]
    for grp in ['T', 'N']:
        g = d[d.grp == grp]; cs = CS.loc[g.index]
        for k in ['stromal', 'immune', 'epithelial']:
            for fp in ['FP:KI_MEK', 'FP:EGF', 'FP:DNA_damage', 'FP:mitosis_cellcycle', 'FP:KI_EGFR', 'FP:KI_MET_ALK', 'FP:KI_CK2', 'FP:metformin', 'PSP:CDK1']:
                if fp in g: q3.append(dict(cohort=coh, grp=grp, covariate=k, fp=fp, rho=spearmanr(g[fp], cs[k], nan_policy='omit')[0]))
    t = CS[(CS.cohort == coh)]
    for k in ['stromal', 'immune', 'epithelial']:
        q3.append(dict(cohort=coh, grp='T-N', covariate=k, fp='level diff (SD units)', rho=t[t.type == 'Tumor'][k].mean() - t[t.type == 'Normal'][k].mean()))
Q3 = pd.DataFrame(q3); Q3.to_csv(TQ + 'T6_composition_correlates.tsv', sep='\t', index=False)
print(Q3[Q3.grp != 'T-N'].groupby(['grp', 'covariate', 'fp']).rho.median().unstack(0).round(2))
print(Q3[Q3.grp == 'T-N'].pivot(index='cohort', columns='covariate', values='rho').round(2))
