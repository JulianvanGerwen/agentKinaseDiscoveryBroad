"""Independent spot-check of KPNA4 S60 in Menyhart 2023 mmc5 (SGC-CK2-1). Untrusted file: read-only via openpyxl, run with python -I."""
import sys, openpyxl
wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
for ws in wb.worksheets:
    hdr = None
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if hdr is None:
            if row and any(isinstance(c, str) and "Gene" in c for c in row if c):
                hdr = [str(c) for c in row]
            continue
        txt = " ".join(str(c) for c in row if c is not None)
        if "KPNA4" in txt and ("60" in txt):
            print(ws.title, [(h, v) for h, v in zip(hdr, row) if v is not None and ("/" in h or "Multiplicity" in h) and "FDR" not in h][:40])
