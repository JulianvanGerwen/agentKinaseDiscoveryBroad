import sys, openpyxl
wb=openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
for ws in wb.worksheets:
    print("### sheet:",ws.title)
    for r in ws.iter_rows(max_row=20000, values_only=True):
        vals=[str(c) for c in r if c is not None and str(c).strip()]
        if vals: print(" | ".join(vals)[:400])
