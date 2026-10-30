"""I10 step 6: annotate ICA phosphomodules.
- cross-study support: module activity (robust z across experiments); n studies with an experiment |z|>=3; top-study share of sum a^2
- kinase enrichment among members vs the ICA site universe: PSP curated kinase (psp_kinases_human), predicted top-1 kinase (motif/ML, kinaseSubMain)
- coarse motif classes from the 15-mer (P+1 proline-directed, R/K at -3 basophilic, acidic +1..+3, pY)
- GO BP (Enrichr 2023 gmt) gene-level hypergeometric; BH within module
Unexplained = no kinase enrichment (PSP or predicted) at q<0.05 and <15% of members with any PSP kinase.
Output: T7_module_annotation.tsv, T7_module_kinase_enrich.tsv, T7_module_go.tsv"""
import sys, numpy as np, pandas as pd, warnings
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
import os; SUF = '_sc' if os.environ.get('ICA_CENTER', '0') == '1' else ''
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
CA = "/workspace/cache/b_I10/"; RT = "/workspace/results/b_I10/tables/"
MOD = pd.read_csv(RT + 'T6_modules' + SUF + '.tsv', sep='\t'); MS = pd.read_csv(RT + 'T6_module_sites' + SUF + '.tsv', sep='\t')
uni = pd.Index(np.load(CA + 'ica_runs' + SUF + '.npz', allow_pickle=True)['sites'])
S = load.sites().set_index('site').reindex(uni)
M = pd.read_csv(CA + 'mod_meta.tsv', sep='\t', dtype={'study': str}).set_index('exp')
act = pd.read_parquet(CA + 'module_activity' + SUF + '.parquet')
w = S.sequence_window.fillna('_' * 15)
motif = pd.DataFrame({'proline': w.str[8] == 'P', 'basophilic': w.str[4].isin(['R', 'K']) | w.str[5].isin(['R', 'K']),
                      'acidophilic': w.str[8:11].apply(lambda s: sum(c in 'DE' for c in s) >= 2), 'pY': S.residue == 'Y'}, index=uni)
def kin_sets(col):
    d = {}
    for s, v in S[col].dropna().items():
        for k in str(v).split(';'):
            if k: d.setdefault(k, set()).add(s)
    return {k: v for k, v in d.items() if len(v) >= 5}
KS = {'psp': kin_sets('psp_kinases_human'), 'pred': kin_sets('pred_top1_kinase')}
go = {}
for line in open('/home/.cache/gseapy/Enrichr.GO_Biological_Process_2023.gmt'):
    p = line.rstrip('\n').split('\t'); go[p[0]] = set(x for x in p[2:] if x)
ugenes = set(S.gene.dropna()); go = {k: v & ugenes for k, v in go.items() if 10 <= len(v & ugenes) <= 300}
rows, krows, grows = [], [], []
for mod in MOD.module:
    mem = MS[MS.module == mod]; ms = set(mem.site); N = len(uni); n = len(ms)
    # activity / cross-study
    a = act[mod].dropna(); az = (a - a.median()) / (1.4826 * (a - a.median()).abs().median())
    st = M.loc[az.index, 'study']
    hit = az.abs() >= 3
    nst = st[hit].nunique(); share = (az ** 2).groupby(st).sum().max() / (az ** 2).sum()
    top = az.abs().sort_values(ascending=False).index[:10]
    best = {}
    for src, sets in KS.items():
        res = []
        for k, s in sets.items():
            x = len(ms & s)
            if x < 3: continue
            o, p = fisher_exact([[x, n - x], [len(s) - x, N - n - len(s) + x]], alternative='greater')
            res.append((k, x, len(s), o, p))
        if res:
            r = pd.DataFrame(res, columns=['kinase', 'k_in_mod', 'k_in_uni', 'OR', 'p']); r['q'] = multipletests(r.p, method='fdr_bh')[1]
            r['module'] = mod; r['source'] = src; krows.append(r)
            best[src] = r.sort_values('p').iloc[0]
    # GO
    genes = set(S.loc[list(ms), 'gene'].dropna()); G = len(ugenes)
    gr = [(t, len(genes & s), len(s), hypergeom.sf(len(genes & s) - 1, G, len(s), len(genes))) for t, s in go.items() if len(genes & s) >= 3]
    gtop = ''
    if gr:
        g = pd.DataFrame(gr, columns=['term', 'k', 'size', 'p']); g['q'] = multipletests(g.p, method='fdr_bh')[1]; g['module'] = mod
        grows.append(g[g.q < 0.1]); g = g.sort_values('p')
        gtop = '; '.join(f"{t.split(' (GO')[0]}({q:.1g})" for t, q in zip(g.term[:3], g.q[:3]))
    mf = motif.loc[list(ms)].mean() / motif.mean()
    fpsp = S.loc[list(ms), 'psp_kinases_human'].notna().mean()
    kq = min([b.q for b in best.values()] + [1])
    rows.append(dict(module=mod, n_sites=n, n_genes=len(genes), n_studies_active=nst, top_study_share=share,
                     top_exps=';'.join(f"{M.loc[e,'treat']}[{M.loc[e,'study']}]{az[e]:+.1f}" for e in top[:6]),
                     psp_best=f"{best['psp'].kinase}(OR {best['psp'].OR:.1f},q {best['psp'].q:.1g})" if 'psp' in best else '',
                     pred_best=f"{best['pred'].kinase}(OR {best['pred'].OR:.1f},q {best['pred'].q:.1g})" if 'pred' in best else '',
                     min_kin_q=kq, frac_psp_annotated=fpsp,
                     motif_enrich=';'.join(f'{k}:{v:.2f}' for k, v in mf.items()), go_top=gtop))
A = pd.DataFrame(rows).merge(MOD[['module', 'n_runs', 'stability']], on='module')
A['cross_study'] = (A.n_studies_active >= 3) & (A.top_study_share < 0.5)
A['unexplained'] = (A.min_kin_q >= 0.05) & (A.frac_psp_annotated < 0.15)
A.to_csv(RT + 'T7_module_annotation' + SUF + '.tsv', sep='\t', index=False)
pd.concat(krows).to_csv(RT + 'T7_module_kinase_enrich' + SUF + '.tsv', sep='\t', index=False)
if grows: pd.concat(grows).to_csv(RT + 'T7_module_go' + SUF + '.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_colwidth', 120)
print(A.drop(columns=['top_exps']).round(3).to_string())
print(A[['module', 'top_exps']].to_string())
print('cross-study', A.cross_study.sum(), 'unexplained', A.unexplained.sum(), 'both', (A.cross_study & A.unexplained).sum())
