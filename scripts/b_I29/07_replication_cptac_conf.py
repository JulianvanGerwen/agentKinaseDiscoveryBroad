"""I29 replication 2: CPTAC LUAD Confirmatory (PDC000490 phosphosite / PDC000489 proteome, CDAP log ratio vs pooled ref).
EGFR/KRAS/NRAS/BRAF/NF1 status from GDC open ssm_occurrences (CPTAC-3). Discovery cases (b_I27 LUAD 110) excluded.
Reads untrusted files as data only (run with python -I)."""
import json, re, sys, numpy as np, pandas as pd, statsmodels.api as sm
DL = '/workspace/cache/b_I29/dl/'; OUT = '/workspace/results/b_I29/'
def z(s): s = pd.to_numeric(s, errors='coerce'); return (s - s.mean()) / s.std()
# aliquot -> case, tumour/normal
cases = json.load(open(DL + 'gdc_lung_aliquots.json'))['data']['hits']
al = {}
for c in cases:
    if c.get('disease_type') != 'Adenomas and Adenocarcinomas': continue
    for s in c.get('samples', []):
        for p in s.get('portions', []) or []:
            for a in p.get('analytes', []) or []:
                for q in a.get('aliquots', []) or []:
                    al[q['submitter_id'][:10]] = (c['submitter_id'], s.get('tissue_type'), s.get('time_between_excision_and_freezing'))
def read_rows(f, genes):
    with open(f) as fh:
        hdr = fh.readline().rstrip('\n').split('\t'); gi = hdr.index('Gene'); rows = {}
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[gi] in genes: rows[p[0]] = p
    D = pd.DataFrame.from_dict(rows, orient='index', columns=hdr)
    return D, hdr
P, hp = read_rows(DL + 'CPTAC3_LUAD_Confirmatory_Phosphoproteome.phosphosite.tmt11.tsv', {'MAPK1', 'MAPK3', 'RPS6KA1', 'RPS6KA3'})
Q, hq = read_rows(DL + 'CPTAC3_LUAD_Confirmatory_Proteome.tmt11.tsv', {'MAPK1', 'MAPK3', 'DUSP6', 'DUSP4', 'EGFR', 'ETV4', 'ETV5', 'SPRED1'})
def mat(D, keycol=None):
    cols = [c for c in D.columns if c.startswith('Log CPT')]
    M = D[cols].replace('', np.nan).astype(float)
    M.columns = [re.search(r'CPT\d{7}', c).group(0)[:10] for c in cols]
    if keycol: M.index = D[keycol]
    return M
PM = mat(P); PM.index = [f"{g}:{i.split(':')[1]}" for g, i in zip(P.Gene, P.index)]
QM = mat(Q, 'Gene')
print(PM.index.tolist())
disc = set(pd.read_parquet('/workspace/results/b_I27/tables/alleles.parquet').query('cohort=="LUAD"').case_id)
meta = pd.DataFrame({k: dict(case=v[0], tissue=v[1], excision=v[2]) for k, v in al.items()}).T
def tum(M):
    keep = [i for i, c in enumerate(M.columns) if c in meta.index and meta.loc[c, 'tissue'] == 'Tumor' and meta.loc[c, 'case'] not in disc]
    T = M.iloc[:, keep].T; T.index = [meta.loc[c, 'case'] for c in T.index]
    T = T.loc[:, ~T.columns.duplicated()]
    return T.groupby(level=0).mean()
PT = tum(PM); QT = tum(QM)
print('confirmatory tumours with phospho:', len(PT), 'protein:', len(QT))
# mutations
h = json.load(open(DL + 'gdc_ssm.json'))['data']['hits']
mut = {}
for x in h:
    cs = x['case']['submitter_id']
    for t in x['ssm']['consequence']:
        tr = t['transcript']
        if tr.get('is_canonical') and tr.get('aa_change'): mut.setdefault(cs, set()).add(tr['gene']['symbol'] + ':' + tr['aa_change'])
A = PT.join(QT.add_prefix('prot_'), how='left')
def status(c):
    m = mut.get(c, set()); g = {x.split(':')[0] for x in m}
    return 'EGFR' if 'EGFR' in g else ('RAS' if g & {'KRAS', 'NRAS', 'BRAF', 'NF1'} else 'WT')
A['group'] = [status(c) for c in A.index]
A['EGFR_alleles'] = [';'.join(sorted(x for x in mut.get(c, set()) if x.startswith('EGFR'))) for c in A.index]
A['excision'] = [pd.to_numeric(meta[meta.case == c].excision, errors='coerce').median() for c in A.index]
print(A.group.value_counts())
for s in ['MAPK3:t202y204', 'MAPK1:t185y187', 'MAPK3:y204', 'MAPK1:y187', 'MAPK3:t202', 'MAPK1:t185']:
    if s in A:
        g = s.split(':')[0]; p = 'prot_' + g
        d = A[[s, p]].dropna()
        if len(d) > 10: b = np.polyfit(d[p], d[s], 1); A[s + '_padj'] = A[s] - (b[0] * A[p] + b[1])
pp = [c for c in ['MAPK3:t202y204_padj', 'MAPK1:t185y187_padj'] if c in A]
A['pERK_pTpY_padj'] = pd.concat([z(A[c]) for c in pp], axis=1).mean(axis=1) if pp else np.nan
A.to_csv(OUT + 'tables/T7_cptac_conf_tumour_table.tsv', sep='\t')
rows = []
for y in [c for c in A.columns if c not in ('group', 'EGFR_alleles', 'excision')]:
    for lab, cov in [('base', []), ('+log excision', ['excision'])]:
        d = A[[y, 'group'] + cov].dropna()
        if d[y].notna().sum() < 15 or (d.group == 'WT').sum() < 4: continue
        X = pd.DataFrame({g: (d.group == g) * 1.0 for g in ['EGFR', 'RAS']}, index=d.index)
        for c in cov: X[c] = z(np.log(d[c].clip(lower=1)))
        Y = (d[y] - A[y].mean()) / A[y].std(); r = sm.OLS(Y, sm.add_constant(X)).fit(cov_type='HC3')
        for g in ['EGFR', 'RAS']:
            rows.append(dict(outcome=y, model=lab, contrast=f'{g}_vs_WT', beta=r.params[g], lo=r.conf_int().loc[g, 0], hi=r.conf_int().loc[g, 1], p=r.pvalues[g], n_g=int((d.group == g).sum()), n_ref=int((d.group == 'WT').sum())))
T = pd.DataFrame(rows); T.to_csv(OUT + 'tables/T8_replication_cptac_conf.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(T.round(3).to_string())
