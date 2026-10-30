import sys, openpyxl, statistics
targets = [
 ("U2SURP","S",37,"HMDASGPSDSDMPSR"),("TRAF4","S",426,"KPGTWRGSLDESSLG"),
 ("PAK2","T",169,"TEAPAVVTEEEDDDE"),("CDV3","S",106,"RVQAMQISSEKEEDD"),
 ("KPNA4","S",60,"HEDICEDSDIDGDYR"),("PHAX","S",368,"ALASLDESQEGHAEA"),
 ("NFXL1","S",331,"ENPCHAGSCQPCPRV"),("ANLN","S",54,"ASNQQPLSGGEEKSC"),
 ("SF3A1","S",451,"RSIREKQSDDEVYAP"),("GPATCH2","S",195,"CRDQDMDSDRAYQYQ"),
 ("CEP152","T",1241,"ELQTLCKTPPRSLSA"),("YTHDC1","S",308,"RKRARGISPIVFDRS"),
 ("FXR1","S",409,"NSELSNPSETESERK"),("ZC3H8","S",77,"SKDYDVYSDNDICSQ"),
]
def f(x):
    try:
        v=float(x)
        return None if v!=v else v
    except: return None

for path,label in [(sys.argv[1],"CX-4945"),(sys.argv[2],"SGC-CK2-1")]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sh in ["All Sites 4h","All Sites 24h"]:
        ws = wb[sh]
        rows=[]
        hdr=None
        for r in ws.iter_rows(values_only=True):
            if hdr is None:
                if r[0]=="Gene Names": hdr=list(r)
                continue
            rows.append(r)
        idx={h:i for i,h in enumerate(hdr) if h}
        print("#"*80); print(f"{label} | {sh} | file={path.split('/')[-1]} | n_rows={len(rows)}")
        # sanity: CK2 motif / PSP CK2 sites
        ml=[]; ml_ck2=[]; ml_psp=[]; hm_ck2=[]
        for r in rows:
            v=f(r[idx['Mean M/L']])
            if v is None: continue
            ml.append(v)
            if str(r[idx['CK2 Motif (HPRD)']])=="True": ml_ck2.append(v)
            if 'CK2' in str(r[idx['PSP Kinase']]).upper():
                ml_psp.append(v)
                h=f(r[idx['Mean H/M']])
                if h is not None: hm_ck2.append(h)
        def s(x): return f"n={len(x)} median={statistics.median(x):+.3f}" if x else "n=0"
        print(f"  SANITY all sites M/L: {s(ml)}")
        print(f"  SANITY CK2-motif(HPRD) M/L: {s(ml_ck2)}")
        print(f"  SANITY PSP CK2 substrate M/L: {s(ml_psp)}  | their H/M (rescue): {s(hm_ck2)}")
        print("  gene\tsite\tmult\tM/L(inh/DMSO,WT)\tH/L(TM+inh/DMSO)\tH/M(rescue)\tsigM/L\tsigH/M\tFDR_M/L\tFDR_H/M\tCK2motif\tPSP\tReg\tRescued\twindow")
        for g,aa,pos,win in targets:
            hits=[r for r in rows if str(r[idx['Gene Names']]).split(';')[0]==g and str(r[idx['Amino Acid']])==aa and str(pos) in str(r[idx['Positions Within Proteins']]).split(';')]
            if not hits:
                hits=[r for r in rows if win[:7] in str(r[idx['Sequence Window']])]
            if not hits:
                print(f"  {g} {aa}{pos}\tNOT QUANTIFIED/NOT PRESENT"); continue
            for r in hits:
                print("  "+"\t".join([g,f"{aa}{pos}",str(r[idx['Multiplicity']]),
                  str(r[idx['Mean M/L']])[:8],str(r[idx['Mean H/L']])[:8],str(r[idx['Mean H/M']])[:8],
                  str(r[idx['T-test Significant M/L']]),str(r[idx['T-test Significant H/M']]),
                  str(r[idx['T-test M/L BH FDR']])[:8],str(r[idx['T-test H/M BH FDR']])[:8],
                  str(r[idx['CK2 Motif (HPRD)']]),str(r[idx['PSP Kinase']]),
                  str(r[idx['Regulation']]),str(r[idx['Rescued']]),str(r[idx['Sequence Window']])]))
    wb.close()
