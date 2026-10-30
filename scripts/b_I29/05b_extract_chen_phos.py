"""I29: stream-extract MAPK1/3 (and RSK, DUSP6, EGFR) phospho rows from Chen 2020 S1F/S1G by UniProt accession. python -I; data only."""
import sys, openpyxl, pandas as pd
src, out = sys.argv[1], sys.argv[2]
ACC = ('P27361', 'P28482', 'Q15418', 'P51812', 'P00533', 'Q16828', 'Q13115', 'P29353', 'Q06124', 'Q13480')
wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
for sh in ['Table S1F_PhosPepLog2TN', 'Table S1G_PhosSiteLog2TN']:
    rows = []
    for i, r in enumerate(wb[sh].iter_rows(values_only=True)):
        txt = ' '.join(str(v) for v in r[:5] if v is not None)
        if i < 3 or any(a in txt for a in ACC): rows.append(r)
    pd.DataFrame(rows).to_csv(f'{out}/{sh.split("_")[0].replace(" ", "")}_acc.tsv', sep='\t', index=False, header=False); print(sh, len(rows), flush=True)
