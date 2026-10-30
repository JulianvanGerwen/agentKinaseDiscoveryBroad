import sys, openpyxl
for p in sys.argv[1:]:
    print("="*70); print(p)
    try:
        wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    except Exception as e:
        print("ERR", e); continue
    for ws in wb.worksheets:
        print(f"  sheet: {ws.title}  dims={ws.calculate_dimension()}")
        for i, row in enumerate(ws.iter_rows(max_row=3, values_only=True)):
            print("   ", [str(c)[:40] if c is not None else '' for c in row[:25]])
            if i>=2: break
    wb.close()
