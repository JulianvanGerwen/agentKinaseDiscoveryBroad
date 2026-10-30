import sys, openpyxl, statistics
targets=[("U2SURP","S",37),("TRAF4","S",426),("PAK2","T",169),("CDV3","S",106),("KPNA4","S",60),
("PHAX","S",368),("NFXL1","S",331),("ANLN","S",54),("SF3A1","S",451),("GPATCH2","S",195),
("CEP152","T",1241),("YTHDC1","S",308),("FXR1","S",409),("ZC3H8","S",77)]
def f(x):
    try:
        v=float(x); return None if v!=v else v
    except: return None
out=open(sys.argv[3],"w")
out.write("dataset\tinhibitor\tdose\ttime\tcell_line\tgene\tsite\tmultiplicity\tlog2_MoverL_inhVsDMSO_WT\tFDR_M_L\tsig_M_L\tlog2_HoverL_TMinh_vs_DMSO\tlog2_HoverM_rescue\tFDR_H_M\tsig_H_M\tCK2motif_HPRD\tRegulation\tRescued_flag\tsequence_window\n")
for path,inh,dose in [(sys.argv[1],"CX-4945","30 uM"),(sys.argv[2],"SGC-CK2-1","5 uM")]:
    wb=openpyxl.load_workbook(path,read_only=True,data_only=True)
    for sh,tm in [("All Sites 4h","4h"),("All Sites 24h","24h")]:
        ws=wb[sh]; hdr=None; rows=[]
        for r in ws.iter_rows(values_only=True):
            if hdr is None:
                if r[0]=="Gene Names": hdr=list(r)
                continue
            rows.append(r)
        I={h:i for i,h in enumerate(hdr) if h}
        down=[(f(r[I['Mean H/M']]),r) for r in rows if (f(r[I['Mean M/L']]) or 0)<-0.585 and str(r[I['T-test Significant M/L']])=='True']
        hm=[d for d,_ in down if d is not None]
        print(f"{inh} {tm}: n_down(sig,M/L<-0.585)={len(down)}  median H/M of down sites={statistics.median(hm):+.3f}" if hm else "")
        for g,aa,pos in targets:
            hits=[r for r in rows if str(r[I['Gene Names']]).split(';')[0]==g and str(r[I['Amino Acid']])==aa and str(pos) in str(r[I['Positions Within Proteins']]).split(';')]
            if not hits:
                out.write(f"Menyhart2023\t{inh}\t{dose}\t{tm}\tU2OS FT CSNK2A1-WT/TM\t{g}\t{aa}{pos}\tNA\tNOT_IN_TABLE\t"+"\t".join(["NA"]*10)+"\n"); continue
            for r in hits:
                ml=f(r[I['Mean M/L']])
                if ml is None:
                    out.write(f"Menyhart2023\t{inh}\t{dose}\t{tm}\tU2OS FT CSNK2A1-WT/TM\t{g}\t{aa}{pos}\t{r[I['Multiplicity']]}\tNOT_QUANTIFIED\t"+"\t".join(["NA"]*8)+f"\t{r[I['Sequence Window']]}\n"); continue
                pct=100*sum(1 for x in hm if x < (f(r[I['Mean H/M']]) or -99))/len(hm) if hm else float('nan')
                out.write("\t".join(["Menyhart2023",inh,dose,tm,"U2OS FT CSNK2A1-WT/TM",g,f"{aa}{pos}",str(r[I['Multiplicity']]),
                  f"{ml:.4f}",str(r[I['T-test M/L BH FDR']])[:10],str(r[I['T-test Significant M/L']]),
                  f"{f(r[I['Mean H/L']]):.4f}" if f(r[I['Mean H/L']]) is not None else "NA",
                  f"{f(r[I['Mean H/M']]):.4f}" if f(r[I['Mean H/M']]) is not None else "NA",
                  str(r[I['T-test H/M BH FDR']])[:10],str(r[I['T-test Significant H/M']]),
                  str(r[I['CK2 Motif (HPRD)']]),str(r[I['Regulation']]),str(r[I['Rescued']]),str(r[I['Sequence Window']])])+"\n")
                if g=="KPNA4": print(f"   KPNA4 S60 H/M percentile among sig-down sites: {pct:.1f}th")
    wb.close()
out.close()
