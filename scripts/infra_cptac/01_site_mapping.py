#!/usr/bin/env python3
"""Map CPTAC BCM-harmonised phosphosite IDs (ENSG|ENSP|S123|15mer|flag; GENCODE v34 positions) to canonical
UniProt site keys `{UniProt}_{Res}{Pos}`.

Strategy
 1. ENSP (versioned) -> SwissProt accession via the BCM 'Gene_annotation_and_representable_isoform_mapping_table'
    (column SwissProt). Fallback: BCM gene_name -> reviewed UniProt entry whose primary gene name matches (unique).
 2. Locate the 15-mer (+/-7 aa, '_' = protein terminus) in the canonical UniProt sequence (reviewed fasta, 'sp|ACC|').
    Exact full window first (choose the hit closest to the GENCODE position if several); then progressively shorter
    centred windows (+/-5, +/-3) but only accepting hits within 30 aa of the GENCODE position (or unique when the
    window is >= 11 aa). The residue letter at the mapped position is checked against the site residue.
 3. Same-position shortcut: if canonical seq[pos-1] == residue and the +/-3 window matches at that position, take it.
Output: /workspace/results/infra_cptac/mapping/site_map.parquet (one row per CPTAC site ID, all cohorts) + stats.
"""
import gzip, re, glob, json, os
import pandas as pd

SRC = '<HACKATHON_ROOT>/kinaseSubMain/external'
V12 = f'{SRC}/cptac_pancan_v1.2'
OUT = '/workspace/results/infra_cptac/mapping'; os.makedirs(OUT, exist_ok=True)
MAPTAB = '/workspace/cache/cptac/bcm_readme_phospho/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt'

# ---- canonical UniProt sequences
seqs = {}
with gzip.open(f'{SRC}/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz', 'rt') as fh:
    acc, buf = None, []
    for line in fh:
        if line.startswith('>'):
            if acc: seqs[acc] = ''.join(buf)
            m = line[1:].split('|'); acc = m[1] if '-' not in m[1] else None; buf = []
        elif acc: buf.append(line.strip())
    if acc: seqs[acc] = ''.join(buf)
print('canonical reviewed sequences', len(seqs))
ann = pd.read_csv(f'{SRC}/uniprot/human_tsv/UP000005640_annotations.tsv.gz', sep='\t', usecols=['Entry', 'Reviewed', 'Gene Names (primary)'])
ann = ann[ann.Reviewed.eq('reviewed')]
g2u = ann.groupby('Gene Names (primary)')['Entry'].agg(list)
g2u = {g: v[0] for g, v in g2u.items() if len(v) == 1}

mt = pd.read_csv(MAPTAB, sep='\t', usecols=['gene', 'protein', 'gene_name', 'SwissProt'], dtype=str)
p2sp = mt.dropna(subset=['protein']).drop_duplicates('protein').set_index('protein')['SwissProt'].to_dict()
ensg2name = mt.drop_duplicates('gene').set_index('gene')['gene_name'].to_dict()

# ---- collect site IDs from all cohorts
ids = set()
for f in glob.glob(f'{V12}/*/*_phospho_site_abundance_*.txt'):
    with open(f) as fh:
        next(fh)
        ids.update(l.split('\t', 1)[0] for l in fh)
ids = sorted(ids)
print('unique CPTAC site IDs', len(ids))


def locate(seq, win, res, pos):
    """Return (pos, method) or (None, reason)."""
    c = 7
    for half, method, need_near in ((7, 'win15', False), (5, 'win11', False), (3, 'win7', True)):
        left = win[c - half:c].lstrip('_'); right = win[c + 1:c + 1 + half].rstrip('_')
        pat = left + res + right
        hits = [m.start() + len(left) + 1 for m in re.finditer(f'(?={re.escape(pat)})', seq)]
        if win[c - half:c].startswith('_') and hits:      # N-terminal padding: must start at protein start
            hits = [h for h in hits if h - len(left) == 1] or hits
        if not hits:
            continue
        best = min(hits, key=lambda h: abs(h - pos))
        if need_near and abs(best - pos) > 30:
            continue
        if len(hits) > 1 and abs(best - pos) > 30:
            continue                                       # ambiguous, far from expected
        return best, method + ('_multi' if len(hits) > 1 else '')
    return None, 'window_not_found'


rows = []
for sid in ids:
    ensg, ensp, rp, win, flag = sid.split('|')
    res, pos = rp[0], int(rp[1:])
    acc = p2sp.get(ensp)
    src = 'ensp_swissprot'
    if not isinstance(acc, str) or acc not in seqs:
        acc = g2u.get(ensg2name.get(ensg)); src = 'gene_name'
    if not isinstance(acc, str) or acc not in seqs:
        rows.append((sid, ensg, ensp, ensg2name.get(ensg), res, pos, win, int(flag), None, None, 'no_uniprot', None)); continue
    seq = seqs[acc]
    if win[7] != res:
        rows.append((sid, ensg, ensp, ensg2name.get(ensg), res, pos, win, int(flag), acc, None, 'window_center_mismatch', src)); continue
    l5, r5 = win[2:7].lstrip('_'), win[8:13].rstrip('_')
    if pos - 1 - len(l5) >= 0 and seq[pos - 1 - len(l5):pos + len(r5)] == l5 + res + r5:
        p, meth = pos, 'same_pos'
    else:
        p, meth = locate(seq, win, res, pos)
    rows.append((sid, ensg, ensp, ensg2name.get(ensg), res, pos, win, int(flag), acc, p, meth, src))

m = pd.DataFrame(rows, columns=['cptac_id', 'ensg', 'ensp', 'gene', 'residue', 'pos_gencode', 'window15', 'isoform_flag',
                                'uniprot', 'pos_uniprot', 'method', 'acc_source'])
ok = m.pos_uniprot.notna()
m['pos_uniprot'] = m.pos_uniprot.astype('Int64')
m['site_key'] = None
m.loc[ok, 'site_key'] = m.loc[ok, 'uniprot'] + '_' + m.loc[ok, 'residue'] + m.loc[ok, 'pos_uniprot'].astype(str)
# validation: residue letter at mapped position in canonical sequence
m['residue_check'] = [bool(k is not None and seqs[u][p - 1] == r) for k, u, p, r in zip(m.site_key, m.uniprot, m.pos_uniprot, m.residue)]
assert m.loc[ok, 'residue_check'].all()
m['pos_shift'] = (m.pos_uniprot - m.pos_gencode)
m.to_parquet(f'{OUT}/site_map.parquet', index=False)
st = {'n_ids': int(len(m)), 'n_mapped': int(ok.sum()), 'rate': float(ok.mean()),
      'n_unique_site_keys': int(m.site_key.nunique()),
      'method': m.method.value_counts().to_dict(), 'acc_source': m[ok].acc_source.value_counts().to_dict(),
      'residue_check_pass_of_mapped': float(m.loc[ok, 'residue_check'].mean()),
      'pos_shift_nonzero': int((m.pos_shift.fillna(0) != 0).sum()),
      'rate_primary_isoform': float(ok[m.isoform_flag == 1].mean()), 'rate_secondary_isoform': float(ok[m.isoform_flag == 2].mean())}
json.dump(st, open(f'{OUT}/site_map_stats.json', 'w'), indent=1)
print(json.dumps(st, indent=1))
