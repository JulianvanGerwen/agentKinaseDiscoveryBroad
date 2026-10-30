"""Light extraction (no re-analysis): PSP-annotated CK2 substrates in Menyhart 2023 tables.
Reproduces the audit_F sanity check (scripts/audit_F/t2_query_menyhart.py) and saves per-site values
so the 'known CK2 substrates' reference in fig3b comes from a saved file.
Usage: python -I s00_menyhart_psp_ck2.py <mmc4.xlsx CX-4945> <mmc5.xlsx SGC-CK2-1> <out.tsv>
"""
import sys, openpyxl, csv
def f(x):
    try:
        v = float(x); return None if v != v else v
    except Exception:
        return None
rows_out = []
for path, label in [(sys.argv[1], "CX-4945"), (sys.argv[2], "SGC-CK2-1")]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sh, t in [("All Sites 4h", "4h"), ("All Sites 24h", "24h")]:
        hdr = None
        for r in wb[sh].iter_rows(values_only=True):
            if hdr is None:
                if r[0] == "Gene Names":
                    hdr = list(r); idx = {h: i for i, h in enumerate(hdr) if h}
                continue
            if 'CK2' not in str(r[idx['PSP Kinase']]).upper():
                continue
            ml = f(r[idx['Mean M/L']])
            if ml is None:
                continue
            rows_out.append([label, t, str(r[idx['Gene Names']]).split(';')[0], str(r[idx['Amino Acid']]),
                             str(r[idx['Positions Within Proteins']]).split(';')[0], str(r[idx['PSP Kinase']]),
                             ml, f(r[idx['Mean H/L']]) if 'Mean H/L' in idx else None, f(r[idx['Mean H/M']])])
with open(sys.argv[3], "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["inhibitor", "time", "gene", "aa", "pos", "psp_kinase", "log2_M_L_WT_inh_vs_DMSO", "log2_H_L_TM_inh_vs_DMSO", "log2_H_M_rescue"])
    w.writerows(rows_out)
print(len(rows_out))
