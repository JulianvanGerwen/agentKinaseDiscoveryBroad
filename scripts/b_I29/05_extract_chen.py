"""I29: stream-extract needed rows from Chen 2020 Table S1 (untrusted xlsx; run with python -I; reads data only).
Usage: python -I 05_extract_chen.py <xlsx> <outdir>"""
import sys, openpyxl, pandas as pd
src, out = sys.argv[1], sys.argv[2]
GENES = {'MAPK1', 'MAPK3', 'ETV4', 'ETV5', 'ETV1', 'SPRED1', 'SPRED2', 'DUSP6', 'DUSP4', 'DUSP5', 'SPRY2', 'SPRY4', 'EGFR', 'RPS6KA1', 'RPS6KA3', 'PHLDA1', 'FOS', 'EGR1', 'KRAS', 'NRAS', 'BRAF', 'HRAS', 'NF1', 'TP53'}
wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
for sh in ['Table S1A_clinical_103patient', 'Table S1C_SNV', 'Table S1D_transcriptome_log2TN', 'Table S1E_ProteomeLog2TN', 'Table S1F_PhosPepLog2TN', 'Table S1G_PhosSiteLog2TN']:
    ws = wb[sh]; rows = []; hdr = None; nread = 0
    for r in ws.iter_rows(values_only=True):
        nread += 1
        if nread <= 4 or sh.startswith('Table S1A'):
            rows.append(r); continue
        txt = ' '.join(str(v) for v in r[:6] if v is not None)
        toks = set(txt.replace(';', ' ').replace('|', ' ').replace('_', ' ').replace(':', ' ').split())
        if toks & GENES: rows.append(r)
    pd.DataFrame(rows).to_csv(f'{out}/{sh.split("_")[0].replace(" ", "")}.tsv', sep='\t', index=False, header=False)
    print(sh, nread, len(rows), flush=True)
