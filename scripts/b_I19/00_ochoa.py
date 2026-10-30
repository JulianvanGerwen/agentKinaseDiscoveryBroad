"""Convert Ochoa 2020 functional score + selected features to parquet (cache)."""
import pandas as pd
D='<HACKATHON_ROOT>/kinaseSubMain/external/ochoa2020/'
fs=pd.read_excel(D+'41587_2019_344_MOESM5_ESM.xlsx',engine='openpyxl')
fs.to_parquet('/workspace/results/b_I19/cache/ochoa_fs.parquet')
cols=['uniprot','position','residue','w0_mya','w3_mya','isInterface','isHotspot','is_disopred','domain_uniprot','PWM_max_mss','netpho_max_KIN','quant_top1','pubmed_counts','paxdb_abundance_log10']
f=pd.read_excel(D+'41587_2019_344_MOESM4_ESM.xlsx',sheet_name='annotated_phosphoproteome',usecols=cols,engine='openpyxl',na_values=['NA'])
f.to_parquet('/workspace/results/b_I19/cache/ochoa_feat.parquet')
print(fs.shape,f.shape)
