"""I27 step 3: exploratory FDR-controlled scan of all outcomes x all allele contrasts with n>=8 per arm."""
import numpy as np, pandas as pd, common as C
from scipy.stats import t as tdist

OUT = '/workspace/results/b_I27/tables/'
A = pd.read_parquet(OUT + 'alleles.parquet')
X = pd.read_parquet(OUT + 'outcomes_tc_padj.parquet')
D = X.merge(A, on=['case_id', 'cohort'], how='left')
YC = [c for c in D.columns if c.startswith(('stage_', 'kin_', 'fp_', 'site_'))]
D['SGKmAKT'] = D.site_NDRG1_SGK - D.site_PRAS40_T246
YC = YC + ['SGKmAKT']


def vec_contrast(d, g, Y, cov):
    """HC3 t-tests of group indicator g for every column of Y (NaNs handled per column)."""
    out = []
    Z = np.hstack([np.ones((len(d), 1)), g[:, None], cov])
    Yv = d[Y].values.astype(float)
    for j, y in enumerate(Y):
        m = np.isfinite(Yv[:, j])
        if m.sum() < 10 or g[m].sum() < 5 or (1 - g[m]).sum() < 5:
            out.append((y, np.nan, np.nan, np.nan, int(g[m].sum()), int((1 - g[m]).sum())))
            continue
        Zm, yv = Z[m], Yv[m, j]
        if np.linalg.matrix_rank(Zm) < Zm.shape[1]:
            Zm = Zm[:, :2]
        b, se, _ = C._ols_hc3(yv, Zm)
        dof = max(1, m.sum() - Zm.shape[1])
        tt = b[1] / se[1] if se[1] > 0 else np.nan
        out.append((y, b[1], se[1], 2 * tdist.sf(abs(tt), dof), int(g[m].sum()), int((1 - g[m]).sum())))
    return out


rows = []
specs = [('KRAS_allele', 'KRAS'), ('PIK3CA_class', 'PIK3CAclass'), ('PIK3CA_allele', 'PIK3CA'),
         ('BRAF_class', 'BRAFclass'), ('EGFR_allele', 'EGFR')]
for coh in C.COHORTS:
    d0 = D[D.cohort == coh]
    cov = d0[['purity', 'prolif']].values.astype(float)
    for gcol, gname in specs:
        vc = d0[gcol].value_counts()
        labs = [l for l in vc.index if vc[l] >= 8 and l not in ('multi', 'other')]
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                a, b = labs[i], labs[j]
                if a == 'WT': a, b = b, a      # orient: mutant vs WT
                m = d0[gcol].isin([a, b])
                dd = d0[m]
                g = (dd[gcol] == a).astype(float).values
                for (y, be, se, p, na, nb) in vec_contrast(dd, g, YC, cov[m.values]):
                    rows.append(dict(cohort=coh, gene=gname, a=a, b=b, y=y, beta=be, se=se, p=p,
                                     n_a=na, n_b=nb))
R = pd.DataFrame(rows).dropna(subset=['p'])
R['lo'] = R.beta - 1.96 * R.se
R['hi'] = R.beta + 1.96 * R.se
R['mde'] = 2.8 * R.se
# BH within the whole exploratory family
o = R.p.values.argsort()
q = np.empty(len(R)); n = len(R)
pv = R.p.values[o]
q[o] = np.minimum.accumulate((pv * n / (np.arange(n) + 1))[::-1])[::-1]
R['q'] = np.clip(q, 0, 1)
R.sort_values('p').to_csv(OUT + 'T4_exploratory.tsv', sep='\t', index=False)

# cross-cohort meta of each (gene, a, b, y) pair seen in >=2 cohorts
metas = []
for k, g in R.groupby(['gene', 'a', 'b', 'y']):
    if g.cohort.nunique() < 2: continue
    m = C.meta(g.to_dict('records'))
    if m.get('k', 0) >= 2:
        metas.append(dict(gene=k[0], a=k[1], b=k[2], y=k[3], **m, cohorts=';'.join(g.cohort)))
M = pd.DataFrame(metas)
if len(M):
    o = M.p_re.values.argsort(); n = len(M); q = np.empty(n)
    pv = M.p_re.values[o]
    q[o] = np.minimum.accumulate((pv * n / (np.arange(n) + 1))[::-1])[::-1]
    M['q_re'] = np.clip(q, 0, 1)
    M.sort_values('p_re').to_csv(OUT + 'T5_exploratory_meta.tsv', sep='\t', index=False)
print(len(R), 'tests; min q', R.q.min(), '; hits q<.1:', (R.q < .1).sum())
print(R.nsmallest(25, 'p')[['cohort', 'gene', 'a', 'b', 'y', 'n_a', 'n_b', 'beta', 'p', 'q']].to_string())
if len(M):
    print(M.nsmallest(15, 'p_re')[['gene', 'a', 'b', 'y', 'k', 'beta_re', 'p_re', 'q_re', 'cohorts']].to_string())
