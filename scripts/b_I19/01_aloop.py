"""Activation-segment (DFG..APE) coordinates for human kinases via the KinCore alignment.
Reference columns taken from AKT1 (DFG D292, APE E319). Output: cache/aloop.tsv (gene, uniprot, dfg_pos, ape_pos)."""
import re, pandas as pd
E='<HACKATHON_ROOT>/kinaseSubMain/external/kincore/'
seqs={};hdr={};name=None
for l in open(E+'alignment-files_Human-PK-alignment.fasta'):
    l=l.strip()
    if l.startswith('>'): name=l[1:].split()[0]; hdr[name]=l[1:].split(); seqs[name]=''
    else: seqs[name]+=l
def colmap(name):
    start=int(name.split('/')[1].split('-')[0]); s=seqs[name]; pos=start-1; m={}
    for i,c in enumerate(s):
        if c.isalpha() and c.isupper(): pos+=1; m[i]=pos
        elif c.isalpha(): pos+=1  # lowercase insert residues still count
    return m
ak=[k for k in seqs if k.startswith('AGC_AKT1/')][0]; m=colmap(ak); inv={v:k for k,v in m.items()}
full=''.join(c for c in seqs[ak] if c.isalpha()).upper()
print('AKT1 290-296', full[290-150:297-150], ' 316-320', full[316-150:321-150])
cD, cE = inv[292], inv[319]
rows=[]
for k,s in seqs.items():
    if k.startswith('0ANN'): continue
    h=hdr[k]; mm=colmap(k)
    def near(c,dirn):
        for d in range(0,15):
            for cc in ([c+d,c-d] if dirn==0 else [c+dirn*d]):
                if cc in mm: return mm[cc]
        return None
    dfg=near(cD,0); ape=near(cE,0)
    if dfg and ape and 15<=ape-dfg<=60:
        rows.append(dict(kc_name=k,gene=h[2] if len(h)>2 else None,uniprot=h[3] if len(h)>3 else None,dfg_pos=dfg,ape_pos=ape))
df=pd.DataFrame(rows); df.to_csv('/workspace/results/b_I19/cache/aloop.tsv',sep='\t',index=False)
print(len(df)); print(df[df.gene.isin(['AKT1','MAP2K1','MAPK1','MAPK3','BRAF','EGFR','MTOR','CDK1','SRC','ABL1','PDPK1','RPS6KB1','GSK3B','MAPK14'])])
