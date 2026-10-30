#!/usr/bin/env python3
"""Download missing CPTAC pan-cancer (Li et al. 2023) open files from PDC (mutation/CNV/RNA) into /workspace/cache/cptac/pdc_dl."""
import json,sys,os,urllib.request,subprocess
out='/workspace/cache/cptac/pdc_dl'; os.makedirs(out,exist_ok=True)
man=json.load(open('<HACKATHON_ROOT>/kinaseSubMain/external/cptac_pancancer_pdc/pancancer_file_manifest.json'))
for name in sys.argv[1:]:
    f=[x for x in man if x['file_name']==name][0]
    u=json.load(urllib.request.urlopen('https://pdc.cancer.gov/pdcapi/file/signedURLFromUuid/'+f['file_id'],timeout=60))['data']
    dest=f"{out}/{name}"
    r=subprocess.run(['curl','-sS','--fail','-L','--retry','3','-o',dest,u]); print(name,r.returncode,os.path.getsize(dest) if os.path.exists(dest) else None,flush=True)
