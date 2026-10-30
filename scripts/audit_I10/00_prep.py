"""audit_I10 step 0: external annotations.
- HPA immune-cell-specific proteins (UniProt): single-cell-type enriched/enhanced/group-enriched where ALL listed specific cell types are
  immune/blood cell types, plus HPA 'RNA blood cell specificity' immune enriched with tissue-cell-type enrichment in immune cells.
- STRING v12 physical links (combined score >= 700) as UniProt pairs (same-complex proxy).
Outputs: cache/audit_I10/{immune_uniprot.txt, string_phys700.parquet}"""
import zipfile, gzip, pandas as pd, numpy as np
EXT = '<HACKATHON_ROOT>/kinaseSubMain/external/'
OUT = '/workspace/cache/audit_I10/'
IMM = ['T-cells', 'B-cells', 'NK-cells', 'monocytes', 'Macrophages', 'Neutrophils', 'cDC', 'pDCs', 'Plasma cells', 'Mast cells',
       'Kupffer cells', 'Microglia', 'Hofbauer cells', 'Innate lymphoid cells', 'Thymocytes', 'Langerhans cells', 'Neutrophil progenitors',
       'Monocyte progenitors', 'dendritic cells', 'Alveolar macrophages', 'Megakaryocytes', 'Platelets', 'Erythroid', 'Basophils', 'Eosinophils',
       'Megakaryocyte progenitors', 'Megakaryocyte-Erythroid progenitors', 'Lymphoid progenitors', 'Granulocytes', 'Erythrocytes']
z = zipfile.ZipFile(EXT + 'hpa/proteinatlas.tsv.zip')
d = pd.read_csv(z.open('proteinatlas.tsv'), sep='\t', usecols=['Gene', 'Uniprot', 'RNA single cell type specificity', 'RNA single cell type specific nCPM',
                                                                  'RNA blood cell specificity', 'RNA tissue cell type enrichment'])
def all_imm(s):
    if not isinstance(s, str): return False
    cts = [x.split(':')[0].strip() for x in s.split(';')]
    return all(any(c.lower().startswith(i.lower()) or i.lower() in c.lower() for i in IMM) for c in cts)
sc = d['RNA single cell type specificity'].isin(['Cell type enriched', 'Group enriched', 'Cell type enhanced']) & d['RNA single cell type specific nCPM'].map(all_imm)
bl = d['RNA blood cell specificity'].isin(['Immune cell enriched', 'Group enriched']) & d['RNA tissue cell type enrichment'].fillna('').str.contains(
    'T-cells|B-cells|macrophages|monocytes|neutrophils|NK-cells|dendritic|plasma|Kupffer|Hofbauer|microglial|Langerhans', case=False)
imm = d[sc | bl].dropna(subset=['Uniprot'])
ids = sorted({u.strip() for s in imm.Uniprot for u in s.split(',')})
open(OUT + 'immune_uniprot.txt', 'w').write('\n'.join(ids))
print('HPA immune-specific proteins', len(ids), 'e.g.', sorted(imm.Gene)[:40])
# STRING physical
al = pd.read_csv(EXT + 'string/9606.protein.aliases.v12.0.txt.gz', sep='\t', names=['sid', 'alias', 'source'], skiprows=1)
al = al[al.source.str.contains('UniProt_AC')]
# keep a single canonical accession per STRING id where possible: use all accessions
m = al.groupby('sid').alias.apply(list)
L = pd.read_csv(EXT + 'string/9606.protein.physical.links.detailed.v12.0.txt.gz', sep=' ', usecols=['protein1', 'protein2', 'combined_score'])
L = L[L.combined_score >= 700]
rows = []
for a, b in zip(L.protein1, L.protein2):
    for x in m.get(a, []):
        for y in m.get(b, []): rows.append((x, y))
P = pd.DataFrame(rows, columns=['u1', 'u2']).drop_duplicates()
P.to_parquet(OUT + 'string_phys700.parquet'); print('STRING physical >=700 uniprot pairs', len(P))
