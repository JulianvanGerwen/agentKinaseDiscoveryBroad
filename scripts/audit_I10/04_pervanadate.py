"""audit_I10 step 4: immune infiltration -> pervanadate (phosphatase_inh) footprint mimicry: leukocyte pY-content confound?
Score = sign-aligned mean of standardised tc_padj sites (own), per cohort OLS: z(score) ~ z(ESTIMATE immune) + purity (median-imputed);
meta = Stouffer z weighted sqrt(n) (as b_I10). Variants of the site set:
  full | no HPA immune-specific proteins | no immune-correlated proteins (protein abundance vs ESTIMATE immune, Spearman >0.3 in >=1 cohort
  or pooled >0.2) | both removed | pS/pT only | pY only | pY only after both removals
Nulls: 200 random site sets matched on coverage x SD x residue (pY vs pS/T) x pro+1, with the footprint's signs (does ANY pY-rich set do this?).
Drivers: per-site meta-z of site ~ immune + purity.
Outputs: A4_pervanadate.tsv, A4_site_drivers.tsv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import statsmodels.api as sm
from scipy.stats import norm, spearmanr
from common import *
rng = np.random.default_rng(9); NN = 200
SA = site_annot(); PV = fp_set('phosphatase_inh')
HPA = set(open(CA + 'immune_uniprot.txt').read().split())
DATA = {}; imm_corr = []
for coh in COH:
    x = cohort(coh); Z = standardise(x)
    ann = cl.table(coh, 'annotations').assign(s=lambda q: q.case_id + '_T').set_index('s').reindex(Z.columns)
    cli = cl.table(coh, 'clinical').assign(s=lambda q: q.case_id + '_T').set_index('s').reindex(Z.columns)
    imm = ann.ESTIMATE_ImmuneScore.astype(float); pur = cli.WES_purity.astype(float); pur = pur.fillna(pur.median())
    pr = cl.matrix(coh, 'protein'); pr = pr[[c for c in Z.columns if c in pr.columns]]; pr = pr[pr.notna().sum(1) >= 20]
    rho = pr.T.corrwith(imm.reindex(pr.columns), method='spearman')
    imm_corr.append(rho.rename(coh))
    DATA[coh] = (Z, x, imm, pur)
IC = pd.concat(imm_corr, axis=1)
EMP = set(IC.index[(IC > 0.3).any(axis=1) | (IC.median(1) > 0.2)])
print('HPA immune', len(HPA), 'empirical immune-correlated proteins', len(EMP), 'overlap', len(HPA & EMP))
prot = lambda idx: np.array([s.split('_')[0] for s in idx]); resid = lambda idx: SA.residue.reindex(idx).fillna('S').values
P_ = prot(PV.index); R_ = resid(PV.index)
VAR = {'full': PV, 'no_HPA_immune': PV[~np.isin(P_, list(HPA))], 'no_immune_corr': PV[~np.isin(P_, list(EMP))],
       'no_HPA_no_immcorr': PV[~np.isin(P_, list(HPA | EMP))], 'pST_only': PV[R_ != 'Y'], 'pY_only': PV[R_ == 'Y'],
       'pY_only_no_HPA_no_immcorr': PV[(R_ == 'Y') & ~np.isin(P_, list(HPA | EMP))]}
def fit(s, imm, pur):
    D = pd.DataFrame({'y': s, 'x': imm, 'p': pur}).dropna()
    if len(D) < 30: return np.nan, np.nan, 0
    D = (D - D.mean()) / D.std()
    r = sm.OLS(D.y, sm.add_constant(D[['x', 'p']])).fit(); return r.params['x'], r.pvalues['x'], len(D)
def meta(per):
    P = pd.DataFrame(per, columns=['b', 'p', 'n']).dropna(); z = norm.isf(P.p / 2) * np.sign(P.b)
    return (z * np.sqrt(P.n)).sum() / np.sqrt(P.n.sum()), int(((P.p < .05) & (P.b > 0)).sum()), P.b.values
rows = []
for vn, sg in VAR.items():
    per, nulls = [], np.zeros((NN, 0)).tolist()
    nullper = [[] for _ in range(NN)]
    for coh, (Z, x, imm, pur) in DATA.items():
        s = score(Z, sg, 8); per.append(fit(s, imm, pur))
        cov = Z.notna().mean(1); sd = x.std(1); rr = resid(Z.index) == 'Y'; p1 = SA.pro1.reindex(Z.index).fillna(False).values
        q = lambda v: pd.qcut(pd.Series(v).rank(method='first'), 4, labels=False).values
        b = ((q(cov.values) * 4 + q(sd.values)) * 2 + rr) * 2 + p1
        pool = {k: np.where(b == k)[0] for k in np.unique(b)}
        sgz = sg[sg.index.isin(Z.index)]; bi = b[[Z.index.get_loc(t) for t in sgz.index]]
        for k in range(NN):
            jj = np.array([rng.choice(pool[v]) for v in bi]); ns = pd.Series(sgz.values, index=Z.index[jj]); ns = ns[~ns.index.duplicated()]
            nullper[k].append(fit(score(Z, ns, 8), imm, pur))
    mz, nrep, bs = meta(per); nz = np.array([meta(npk)[0] for npk in nullper])
    rows.append(dict(variant=vn, n_sites=len(sg), frac_pY=(resid(sg.index) == 'Y').mean(), meta_z=mz, n_pos_sig=nrep,
                     cohort_betas=';'.join(f'{c}:{b:+.2f}' for c, b in zip(DATA, bs)), null_meta_z_median=np.median(nz),
                     null_q95=np.quantile(nz, .95), p_vs_matched_null=(1 + (nz >= mz).sum()) / (NN + 1)))
    print(rows[-1], flush=True)
T = pd.DataFrame(rows); T.to_csv(RA + 'A4_pervanadate.tsv', sep='\t', index=False)
# site drivers
dr = []
for t in PV.index:
    per = []
    for coh, (Z, x, imm, pur) in DATA.items():
        if t in Z.index: per.append(fit(Z.loc[t] * PV[t], imm, pur))
    if len(per) >= 3:
        mz, nrep, bs = meta(per); u = t.split('_')[0]
        dr.append(dict(site=t, gene=SA.gene.get(t, ''), residue=t.split('_')[1][0], sign=PV[t], n_coh=len(per), meta_z=mz, n_pos_sig=nrep,
                       HPA_immune=u in HPA, immune_corr_protein=u in EMP, prot_imm_rho_median=IC.loc[u].median() if u in IC.index else np.nan))
D = pd.DataFrame(dr).sort_values('meta_z', ascending=False); D.to_csv(RA + 'A4_site_drivers.tsv', sep='\t', index=False)
print(D.head(30).round(2).to_string()); print(D.groupby(['HPA_immune', 'immune_corr_protein']).meta_z.agg(['size', 'median', 'mean']).round(2))
print(D.groupby('residue').meta_z.agg(['size', 'median', 'mean']).round(2))
