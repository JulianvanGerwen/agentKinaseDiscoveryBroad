"""Task 3d: what predicts cross-study invariance of responding sites?
Unit: (group, site) with >=3 measured studies and >=1 responding study (up run or down run).
Outcome: replicated = responds in >=2 studies.  Predictors: Ochoa 2020 functional score (NB: partly built from regulation features -> circular-ish),
site age (w0_mya), disorder, PSP regulatory, residue Y, protein is a kinase, kinase-motif class of top predicted kinase,
covariates: n measured studies, group fixed effects.  Logistic regression with site-clustered SE (statsmodels GLM)."""
import sys, glob, numpy as np, pandas as pd
import statsmodels.formula.api as smf, statsmodels.api as sm
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
OUT = "/workspace/results/b_I19/"; CA = OUT + 'cache/'
S = load.sites().set_index('site'); K = load.kinases()
fs = pd.read_parquet(CA + 'ochoa_fs.parquet'); fs['site'] = fs.uniprot + '_' + fs.position.astype(int).astype(str)
ft = pd.read_parquet(CA + 'ochoa_feat.parquet'); ft['site'] = ft.uniprot + '_' + ft.residue + ft.position.astype(int).astype(str)
fs = fs.merge(ft[['site', 'uniprot', 'position']].assign(site2=ft.site), on=['uniprot', 'position'], how='left')
fsd = fs.drop_duplicates('site2').set_index('site2').functional_score
ftd = ft.drop_duplicates('site').set_index('site')
kin_genes = set(K.kinase)
fam = dict(zip(K.kinase, K.manning_group))
def motif_class(k):
    g = fam.get(k, None)
    if k in ('CSNK2A1', 'CSNK2A2', 'PLK1', 'CSNK1A1', 'CSNK1D', 'CSNK1E', 'GSK3B', 'GSK3A'): return 'acidophilic/primed'
    if g in ('CMGC',): return 'proline-directed'
    if g in ('AGC', 'CAMK'): return 'basophilic'
    if g in ('TK',): return 'tyrosine'
    return 'other' if isinstance(k, str) else 'none'
rows = []
for d in ['up', 'down']:
    suf = '' if d == 'up' else '_down'
    for f in glob.glob(CA + f'site_stats{suf}_*.parquet'):
        g = f.split(f'site_stats{suf}_')[1].replace('.parquet', '')
        if suf == '' and g.startswith('down_'): continue
        t = pd.read_parquet(f); t = t[t.nup >= 1].copy(); t['group'] = g; t['dir'] = d; rows.append(t)
T = pd.concat(rows)
T = T[~T.group.str.startswith('AXIS')]  # avoid double-counting with class groups
T['replicated'] = (T.nup >= 2).astype(int)
T['fs'] = T.site.map(fsd); T['age'] = T.site.map(ftd.w0_mya); T['disorder'] = T.site.map(ftd.is_disopred)
T['psp_reg'] = T.site.map(S.psp_regulatory).fillna(False).astype(int)
T['resY'] = (T.site.map(S.residue) == 'Y').astype(int)
T['is_kinase'] = T.site.map(S.gene).isin(kin_genes).astype(int)
T['motif'] = T.site.map(S.pred_top1_kinase).map(motif_class)
print(T.groupby('dir').agg(n=('site', 'size'), rep=('replicated', 'mean'), fs_cov=('fs', lambda x: x.notna().mean()), age_cov=('age', lambda x: x.notna().mean())))
T.to_csv(OUT + 'tables/T3_predictor_table.tsv.gz', sep='\t', index=False)
res = []
for d in ['up', 'down']:
    D = T[(T.dir == d)].dropna(subset=['fs', 'age', 'disorder']).copy()
    D['fs_z'] = (D.fs - D.fs.mean()) / D.fs.std(); D['age_z'] = (D.age - D.age.mean()) / D.age.std()
    for form in ['replicated ~ fs_z + meas + C(group)', 'replicated ~ age_z + disorder + psp_reg + is_kinase + C(motif, Treatment("acidophilic/primed")) + meas + C(group)']:
        m = smf.glm(form, D, family=sm.families.Binomial()).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(D.site)[0]})
        for v in m.params.index:
            if v.startswith('C(group)') or v == 'Intercept': continue
            res.append(dict(dir=d, model=form.split('~')[1].strip()[:40], term=v, OR=np.exp(m.params[v]), lo=np.exp(m.conf_int().loc[v, 0]), hi=np.exp(m.conf_int().loc[v, 1]), p=m.pvalues[v], n=len(D), n_rep=int(D.replicated.sum())))
Rz = pd.DataFrame(res); Rz.to_csv(OUT + 'tables/T3_predictors_logit.tsv', sep='\t', index=False)
pd.set_option('display.width', 200); print(Rz.round(4).to_string())
