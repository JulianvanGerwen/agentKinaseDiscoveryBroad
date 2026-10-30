import os,sys,time
os.environ['OMP_NUM_THREADS']='4'
import pandas as pd, pyarrow.parquet as pq
from threadpoolctl import threadpool_limits
sys.path.insert(0,'/workspace/scripts/infra_compendium'); sys.path.insert(0,'/workspace/scripts/infra_compendium/kinact_mvp_pkg')
import load
from kinact_mvp import KinaseActivityModel
model=KinaseActivityModel.from_config('/workspace/scripts/infra_compendium/kinact_mvp_pkg/model/config.json')
prots=['P00519','Q14004','Q05397','Q13177','O14920','P28482','P27361']
X=load.sites_matrix(); mask=X.index.str.split('_').str[0].isin(prots); X.loc[mask]=float('nan')
kin=['ABL1','CDK13','PTK2','PAK2','IKBKB','MAPK1','MAPK3']
cols=list(X.columns); out=[]; t0=time.time()
for i in range(0,len(cols),250):
    with threadpool_limits(limits=4):
        r=model.score(X[cols[i:i+250]],device='cpu',n_jobs=4)
    lg=r.long(); out.append(lg[lg.kinase.isin(kin)]); print(i,time.time()-t0,flush=True)
pd.concat(out).to_parquet('/workspace/results/val_1/tables/corev3_masked_long.parquet',index=False)
