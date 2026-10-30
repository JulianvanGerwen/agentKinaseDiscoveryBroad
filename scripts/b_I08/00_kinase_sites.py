#!/usr/bin/env python3
"""I08 step 0: per-kinase site annotation. PSP regulatory sites (human, phospho) on kinases with function
'activity, induced' (activating) / 'activity, inhibited' (inhibitory); activation loop (DFG..APE heuristic on
canonical UniProt sequence). Output: results/b_I08/tables/kinase_site_annotation.tsv (one row per measured kinase site)."""
import gzip, re, sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
R = cl.ROOT
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


def aloop(seq):
    """activation segment: from DFG (or DLG/DFA/DYG) +3 to the APE-like motif (xPE / xxE within 15-45 aa)."""
    out = []
    for m in re.finditer(r'D[FLYW]G', seq):
        s = m.start()
        win = seq[s + 15:s + 50]
        mm = re.search(r'[AS]PE|[ATSV]PE|P[EI]', win)
        if mm:
            out.append((s + 4, s + 15 + mm.start() + 1))   # 1-based inclusive range approx
    return out


rows = []
for _, k in kin.iterrows():
    sq = seqs.get(k.uniprot)
    loops = aloop(sq) if sq else []
    # keep the first loop after position 100 typical of a kinase domain (most kinases have one DFG)
    rows.append((k.kinase, k.uniprot, loops[0] if loops else None, len(loops)))
loopdf = pd.DataFrame(rows, columns=['kinase', 'uniprot', 'aloop', 'n_dfg'])

reg = pd.read_csv('/workspace/data/Regulatory_sites', sep='\t', skiprows=3, dtype=str)
reg = reg[(reg.ORGANISM == 'human') & reg.MOD_RSD.str.endswith('-p', na=False)].copy()
reg['site_key'] = reg.ACC_ID + '_' + reg.MOD_RSD.str[:-2]
f = reg.ON_FUNCTION.fillna('')
reg['activating'] = f.str.contains('activity, induced')
reg['inhibitory'] = f.str.contains('activity, inhibited')
reg = reg.groupby('site_key').agg(psp_activating=('activating', 'any'), psp_inhibitory=('inhibitory', 'any'),
                                  psp_function=('ON_FUNCTION', lambda x: '; '.join(sorted(set(x.dropna()))))).reset_index()

# all measured sites on kinases across cohorts
meas = {}
for c in cl.COHORTS:
    idx = cl.matrix(c, 'input_tc').index
    for s in idx:
        a = s.split('_')[0]
        meas.setdefault(s, []).append(c)
acc2k = dict(zip(kin.uniprot, kin.kinase))
sites = pd.DataFrame([(s, s.split('_')[0], int(s.split('_')[1][1:]), s.split('_')[1][0], len(v), ','.join(v))
                      for s, v in meas.items() if s.split('_')[0] in acc2k],
                     columns=['site_key', 'uniprot', 'pos', 'res', 'n_cohorts', 'cohorts'])
sites['kinase'] = sites.uniprot.map(acc2k)
sites = sites.merge(reg, on='site_key', how='left')
for col in ['psp_activating', 'psp_inhibitory']: sites[col] = sites[col].fillna(False).astype(bool)
sites['psp_regulatory'] = sites.psp_function.notna()
lp = dict(zip(loopdf.kinase, loopdf.aloop))
sites['in_aloop'] = [bool(lp.get(k)) and lp[k][0] - 3 <= p <= lp[k][1] for k, p in zip(sites.kinase, sites.pos)]
sites.to_csv('/workspace/results/b_I08/tables/kinase_site_annotation.tsv', sep='\t', index=False)
loopdf.to_csv('/workspace/results/b_I08/tables/kinase_aloop.tsv', sep='\t', index=False)
print(sites.shape, sites.psp_activating.sum(), sites.in_aloop.sum(), (sites.in_aloop & sites.psp_activating).sum())
chk = sites[sites.kinase.isin(['MAPK1', 'MAPK3', 'CDK1', 'CDK2', 'AKT1', 'SRC', 'PRKACA', 'AURKA', 'EGFR'])]
print(chk[chk.in_aloop | chk.psp_activating][['site_key', 'kinase', 'in_aloop', 'psp_activating', 'n_cohorts']].to_string())
