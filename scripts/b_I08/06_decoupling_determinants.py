#!/usr/bin/env python3
"""I08 step 6: what kinase properties predict abundance-coupling vs site-coupling?
Features: Manning group; RD vs non-RD (catalytic-loop HRD motif from sequence); protein variability across tumours
(median SD of log2 protein); #substrates; metabolic/atypical. Tests: Mann-Whitney / Spearman on meta partial r.
Output: tables/decoupling_determinants.tsv, tables/kinase_properties.tsv"""
import os, sys, re, gzip, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from i08util import *
from scipy import stats
kin = pd.read_csv(f'{R}/model_corev3/kinases.tsv', sep='\t')
SRC = '<HACKATHON_ROOT>/kinaseSubMain/external'
seqs = {}
with gzip.open(f'{SRC}/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz', 'rt') as fh:
    acc, buf = None, []
    for line in fh:
        if line.startswith('>'):
            if acc: seqs[acc] = ''.join(buf)
            m = line[1:].split('|'); acc = m[1] if '-' not in m[1] else None; buf = []
        else: buf.append(line.strip())
    if acc: seqs[acc] = ''.join(buf)


def rd(seq):
    # catalytic loop: H-x-D-[LIVM]-[KRAHS]-x-x-N  (eukaryotic protein kinase HRDxKxxN)
    m = re.search(r'[HY](.)D[LIVMF][KRAHSQ].{2}N', seq or '')
    return (m.group(1) == 'R') if m else None


kin['RD'] = [rd(seqs.get(a)) for a in kin.uniprot]
# protein variability
sd = {}
for c in cl.COHORTS:
    P = cl.matrix(c, 'protein')[cl.tumour_ids(c)]
    for a in kin.uniprot:
        if a in P.index and P.loc[a].notna().sum() >= 20: sd.setdefault(a, []).append(P.loc[a].std())
kin['prot_sd'] = kin.uniprot.map(lambda a: np.median(sd[a]) if a in sd else np.nan)
n = cl.activity('tc', 'corev3_n')
kin['n_sub'] = kin.kinase.map(n.median(1))
kin.to_csv(f'{OUT}/tables/kinase_properties.tsv', sep='\t', index=False)
C = pd.read_csv(f'{OUT}/tables/kinase_classification_all_sources.tsv', sep='\t').merge(kin[['kinase', 'RD', 'prot_sd', 'n_sub']], on='kinase')
rows = []
for src in ['corev3_tc_padj', 'meanz_tc_padj', 'corev3_tc', 'meanz_tc']:
    x = C[C.source == src]
    for y in ['prot_rp', 'site_rp']:
        a = x[x.RD == True][y].dropna(); b = x[x.RD == False][y].dropna()
        if len(a) > 3 and len(b) > 3:
            rows.append(dict(source=src, y=y, test='RD vs nonRD', n1=len(a), n2=len(b), med1=a.median(), med2=b.median(), p=stats.mannwhitneyu(a, b).pvalue))
        a = x[x.group == 'Atypical'][y].dropna(); b = x[x.group != 'Atypical'][y].dropna()
        if len(a) > 3: rows.append(dict(source=src, y=y, test='Atypical vs rest', n1=len(a), n2=len(b), med1=a.median(), med2=b.median(), p=stats.mannwhitneyu(a, b).pvalue))
        a = x[x.group == 'TK'][y].dropna(); b = x[x.group != 'TK'][y].dropna()
        if len(a) > 3: rows.append(dict(source=src, y=y, test='TK vs rest', n1=len(a), n2=len(b), med1=a.median(), med2=b.median(), p=stats.mannwhitneyu(a, b).pvalue))
        for f in ['prot_sd', 'n_sub']:
            m = x[[y, f]].dropna()
            if len(m) > 10:
                rr, p = stats.spearmanr(m[y], m[f]); rows.append(dict(source=src, y=y, test=f'spearman vs {f}', n1=len(m), med1=rr, p=p))
    # paired: site vs protein within kinase
    m = x[['prot_rp', 'site_rp']].dropna()
    if len(m) > 5:
        rows.append(dict(source=src, y='site_rp - prot_rp', test='Wilcoxon paired', n1=len(m), med1=(m.site_rp - m.prot_rp).median(),
                         med2=(m.site_rp > m.prot_rp).mean(), p=stats.wilcoxon(m.site_rp, m.prot_rp).pvalue))
D = pd.DataFrame(rows); D.to_csv(f'{OUT}/tables/decoupling_determinants.tsv', sep='\t', index=False)
pd.set_option('display.width', 200); print(D.round(4).to_string())
print(kin.RD.value_counts(dropna=False))
