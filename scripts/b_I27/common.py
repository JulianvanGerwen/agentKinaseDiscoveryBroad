"""I27 common: allele calls, outcome matrix, regression + permutation + meta."""
import sys, re, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_cptac')
import cptac_load as cl

OUT = '/workspace/results/b_I27'
COHORTS = cl.COHORTS
RNG = np.random.default_rng(7)

HELICAL = {'E542K', 'E545K', 'E545A', 'E545G', 'E545Q', 'Q546R', 'Q546P', 'Q546K', 'Q546E'}
KINASE = {'H1047R', 'H1047L', 'H1047Y', 'M1043V', 'M1043I', 'G1049R', 'H1048R', 'N1044S', 'N1044K'}
BRAF1 = {'V600E', 'V600K', 'V600D', 'V600R'}
BRAF2 = {'G469A', 'G469R', 'G469V', 'K601E', 'K601N', 'L597R', 'L597Q', 'L597S', 'V487_P492delinsA'}
BRAF3 = {'D594G', 'D594N', 'D594A', 'G466V', 'G466E', 'G466A', 'N581S', 'N581I', 'G596R', 'G596D'}


def _clean(pc):
    return pc[2:] if isinstance(pc, str) and pc.startswith('p.') else pc


def allele_table():
    """case_id x {gene}_allele labels, plus cohort."""
    rows = []
    tum = cl.samples()
    tum = tum[(tum.sample_type == 'Tumor') & tum.has_phospho]
    for coh in COHORTS:
        cases = sorted(tum[tum.cohort == coh].case_id.unique())
        m = cl.table(coh, 'mutations')
        m = m[m.case_id.isin(cases)].copy()
        m['pc'] = m.protein_change.map(_clean)
        nonsyn = m[m.Variant_Classification != 'Silent']
        for case in cases:
            d = {'case_id': case, 'cohort': coh}
            g = nonsyn[nonsyn.case_id == case]
            for gene in ['KRAS', 'NRAS', 'PIK3CA', 'BRAF', 'EGFR', 'TP53', 'PTEN', 'STK11', 'KEAP1', 'NF1', 'APC']:
                d[gene + '_mut'] = int((g.gene == gene).any())
            # KRAS
            k = sorted(set(g[g.gene == 'KRAS'].pc.dropna()))
            d['KRAS_allele'] = 'WT' if not k else (k[0] if len(k) == 1 else 'multi')
            # PIK3CA
            p = sorted(set(g[(g.gene == 'PIK3CA')].pc.dropna()))
            pm = [x for x in p if re.fullmatch(r'[A-Z]\d+[A-Z]', str(x))]
            if not p:
                cls, al = 'WT', 'WT'
            elif len(pm) >= 2:
                cls, al = 'double', 'double'
            elif len(pm) == 1:
                al = pm[0]
                cls = 'helical' if al in HELICAL else ('kinase' if al in KINASE else 'other')
            else:
                cls, al = 'other', (p[0] if p else 'other')
            d['PIK3CA_allele'], d['PIK3CA_class'] = al, cls
            # BRAF
            b = sorted(set(g[g.gene == 'BRAF'].pc.dropna()))
            if not b:
                d['BRAF_allele'], d['BRAF_class'] = 'WT', 'WT'
            else:
                d['BRAF_allele'] = b[0] if len(b) == 1 else 'multi'
                cls = 'other'
                for x in b:
                    if x in BRAF1: cls = 'class1'
                    elif x in BRAF2 and cls != 'class1': cls = 'class2'
                    elif x in BRAF3 and cls not in ('class1', 'class2'): cls = 'class3'
                d['BRAF_class'] = cls
            # EGFR
            e = sorted(set(g[g.gene == 'EGFR'].pc.dropna()))
            lab = 'WT'
            if e:
                lab = 'other'
                for x in e:
                    mm = re.match(r'[A-Z](\d+)_[A-Z](\d+)del', str(x))
                    if mm and 740 <= int(mm.group(1)) <= 755:
                        lab = 'ex19del'
                    elif x == 'L858R':
                        lab = 'L858R'
                if len(e) > 1 and lab in ('ex19del', 'L858R'):
                    lab = 'multi'
            d['EGFR_allele'] = lab
            rows.append(d)
    return pd.DataFrame(rows)


def outcomes(variant='tc_padj'):
    """sample-level outcome matrix (z within cohort) + covariates, indexed by sample_id."""
    st = pd.read_parquet(f'/workspace/results/b_I24/tables/stages_{variant}.parquet')
    stage_cols = ['RTK', 'EGFRpY', 'ADAPT', 'ERK', 'RSK', 'ERKSUB', 'AKTSUB', 'MTORC1', 'S6']
    X = st[stage_cols + ['cohort', 'case_id', 'purity', 'prolif']].copy()
    X.columns = ['stage_' + c for c in stage_cols] + ['cohort', 'case_id', 'purity', 'prolif']
    # kinase activity
    z = cl.activity(variant, 'corev3_z')
    n = cl.activity(variant, 'corev3_n')
    keep = n.index[(n >= 5).mean(1) > 0.5]
    z = z.loc[keep]
    z = z.loc[z.notna().mean(1) > 0.5]
    K = z.T
    K.columns = ['kin_' + c for c in K.columns]
    X = X.join(K)
    # footprints
    fp = pd.read_parquet(f'/workspace/cache/b_I10/cptac_scores_{variant}.parquet')
    fp = fp[~fp.footprint.str.contains(r'\|')]
    F = fp.pivot(index='sample', columns='footprint', values='nz')
    F.columns = ['fp_' + c for c in F.columns]
    X = X.join(F)
    # pre-specified sites
    sites = {'site_PRAS40_T246': ['Q96B36_T246'], 'site_NDRG1_SGK': ['Q92597_T328', 'Q92597_T366'],
             'site_GSK3B_S9': ['P49841_S9'], 'site_S6_S235': ['P62753_S235', 'P62753_S236']}
    addl = {k: [] for k in sites}
    idx = []
    for coh in COHORTS:
        ph = cl.matrix(coh, f'input_{variant}')
        tt = [s for s in cl.tumour_ids(coh) if s in ph.columns]
        idx += tt
        for k, ss in sites.items():
            ss2 = [s for s in ss if s in ph.index]
            addl[k].append(ph.loc[ss2, tt].mean(0) if ss2 else pd.Series(np.nan, index=tt))
    S = pd.DataFrame({k: pd.concat(v) for k, v in addl.items()})
    X = X.join(S)
    val_cols = [c for c in X.columns if c.startswith(('stage_', 'kin_', 'fp_', 'site_'))]
    X[val_cols] = X.groupby('cohort')[val_cols].transform(lambda s: (s - s.mean()) / s.std(ddof=1))
    for c in ['purity', 'prolif']:
        X[c] = X.groupby('cohort')[c].transform(lambda s: (s - s.mean()) / s.std(ddof=1))
    return X, val_cols


def _ols_hc3(y, Xd):
    XtXi = np.linalg.pinv(Xd.T @ Xd)
    b = XtXi @ Xd.T @ y
    r = y - Xd @ b
    h = np.einsum('ij,jk,ik->i', Xd, XtXi, Xd)
    h = np.clip(h, 0, 0.999)
    S = (Xd * (r / (1 - h))[:, None])
    V = XtXi @ (S.T @ S) @ XtXi
    return b, np.sqrt(np.maximum(np.diag(V), 0)), r


def contrast(df, ycol, gcol, a, b, covars=('purity', 'prolif'), extra=(), nperm=5000, seed=0):
    """effect of group a relative to group b on ycol. Returns dict."""
    d = df[df[gcol].isin([a, b])].copy()
    cv = [c for c in list(covars) + list(extra) if c in d.columns]
    d = d[[ycol, gcol] + cv].dropna()
    na, nb = int((d[gcol] == a).sum()), int((d[gcol] == b).sum())
    base = dict(y=ycol, group=gcol, a=a, b=b, n_a=na, n_b=nb)
    if na < 3 or nb < 3:
        return {**base, 'beta': np.nan, 'se': np.nan, 'lo': np.nan, 'hi': np.nan, 'p': np.nan, 'mde': np.nan}
    g = (d[gcol] == a).astype(float).values
    C = []
    for c in cv:
        v = d[c]
        if not pd.api.types.is_numeric_dtype(v):
            du = pd.get_dummies(v, drop_first=True).astype(float)
            du = du.loc[:, du.sum() > 0]
            if du.shape[1] and du.shape[1] < len(d) - 4:
                C.append(du.values)
        else:
            vv = v.astype(float).values
            if np.nanstd(vv) > 0:
                C.append(vv[:, None])
    C = np.hstack(C) if C else np.zeros((len(d), 0))
    Xd = np.hstack([np.ones((len(d), 1)), g[:, None], C])
    y = d[ycol].astype(float).values
    if np.linalg.matrix_rank(Xd) < Xd.shape[1]:
        Xd = np.hstack([np.ones((len(d), 1)), g[:, None]])
    bhat, se, _ = _ols_hc3(y, Xd)
    beta, sb = bhat[1], se[1]
    # Freedman-Lane: residualise y on covariates, permute residuals
    Z = np.hstack([np.ones((len(d), 1)), C])
    bz, _, rz = _ols_hc3(y, Z)
    fit = y - rz
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(nperm):
        yp = fit + rng.permutation(rz)
        bp, _, _ = _ols_hc3(yp, Xd)
        cnt += abs(bp[1]) >= abs(beta) - 1e-12
    p = (cnt + 1) / (nperm + 1)
    return {**base, 'beta': beta, 'se': sb, 'lo': beta - 1.96 * sb, 'hi': beta + 1.96 * sb,
            'p': p, 'mde': 2.8 * sb}


def meta(rows):
    """rows: list of dicts with beta, se. Returns FE and RE."""
    r = pd.DataFrame(rows).dropna(subset=['beta', 'se'])
    r = r[r.se > 0]
    if not len(r):
        return dict(k=0)
    w = 1 / r.se ** 2
    bfe = (w * r.beta).sum() / w.sum()
    sfe = np.sqrt(1 / w.sum())
    Q = (w * (r.beta - bfe) ** 2).sum()
    k = len(r)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    w2 = 1 / (r.se ** 2 + tau2)
    bre = (w2 * r.beta).sum() / w2.sum()
    sre = np.sqrt(1 / w2.sum())
    from scipy.stats import norm
    return dict(k=k, beta_fe=bfe, se_fe=sfe, p_fe=2 * norm.sf(abs(bfe / sfe)),
                beta_re=bre, se_re=sre, p_re=2 * norm.sf(abs(bre / sre)),
                lo_re=bre - 1.96 * sre, hi_re=bre + 1.96 * sre, tau2=tau2, Q=Q,
                mde_re=2.8 * sre, n_a=int(r.n_a.sum()), n_b=int(r.n_b.sum()))
