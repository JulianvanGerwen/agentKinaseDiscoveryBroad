"""Parse GDSC1/2 + PRISM secondary into one long table (run with python -I). Out: cache/b_I23/sens_long.parquet"""
import pandas as pd, re, numpy as np
DL="/workspace/cache/b_I23/dl/"; OUT="/workspace/cache/b_I23/"
norm=lambda s: re.sub(r'[^A-Z0-9]','',str(s).upper())
rows=[]
for ds in ["GDSC1","GDSC2"]:
    g=pd.read_excel(DL+f"{ds}_fitted_dose_response_27Oct23.xlsx")
    g=g.rename(columns=str.upper)
    rows.append(pd.DataFrame({"source":ds,"drug":g.DRUG_NAME,"drug_id":g.DRUG_ID.astype(str),"cell":g.CELL_LINE_NAME,"cell_norm":g.CELL_LINE_NAME.map(norm),
        "depmap_id":None,"cosmic":g.COSMIC_ID,"auc":g.AUC,"ln_ic50":g.LN_IC50,"max_conc_uM":g.MAX_CONC,"tissue":g.get("TCGA_DESC")}))
p=pd.read_csv(DL+"secondary-screen-dose-response-curve-parameters.csv",low_memory=False)
print(p.columns.tolist())
rows.append(pd.DataFrame({"source":"PRISM","drug":p["name"],"drug_id":p.broad_id.astype(str)+"|"+p.screen_id.astype(str),"cell":p.ccle_name.astype(str).str.split("_").str[0],
    "cell_norm":p.ccle_name.astype(str).str.split("_").str[0].map(norm),"depmap_id":p.depmap_id,"cosmic":None,"auc":p.auc,"ln_ic50":np.log(p.ic50.where(p.ic50>0)),
    "max_conc_uM":None,"tissue":p.ccle_name.astype(str).str.split("_",n=1).str[1]}))
L=pd.concat(rows,ignore_index=True); L["drug_norm"]=L.drug.map(norm)
L["ln_ic50"]=pd.to_numeric(L.ln_ic50,errors="coerce")
L.to_parquet(OUT+"sens_long.parquet"); print(L.groupby("source").size())
sc=pd.read_csv(DL+"screened_compounds_rel_8.5.csv"); sc.to_csv(OUT+"gdsc_compounds.tsv",sep="\t",index=False); print(sc.columns.tolist())
