"""audit_D shared helpers: robust per-condition z for the whole compendium (cached), condition table, drug classes."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load
CA = "/workspace/cache/audit_D"

def zmat():
    f = f"{CA}/z.parquet"
    if os.path.exists(f):
        return pd.read_parquet(f)
    X = load.sites_matrix()
    V = X.values.astype(np.float32)
    med = np.nanmedian(V, 0); mad = np.nanmedian(np.abs(V - med), 0) * 1.4826
    Z = pd.DataFrame(((V - med) / mad).astype(np.float32), index=X.index, columns=X.columns)
    Z.to_parquet(f)
    return Z

def lfc(sites=None, conds=None):
    return load.sites_matrix(conditions=conds, sites=sites)

def conds():
    C = load.conditions()
    C["usable"] = ~C.qc_flag.isin(["duplicate_vector", "scale_suspect"])
    single = (~C.is_combination.astype(bool)) & C.agents.notna() & ~C.agents.fillna("").str.contains(";")
    C["single"] = single
    return C

CLS = {
 "cat_mTORi": "AZD8055 KU0063794 torin1 torin2 INK128 sapanisertib PP242 torkinib mTORi",
 "dual_PI3K_mTORi": "dactolisib omipalisib PI103",
 "rapalog": "rapamycin everolimus temsirolimus",
 "S6Ki": "PF4708671 LY2584702",
 "AKTi": "MK2206 capivasertib ipatasertib GSK690693 akti_akt12 aktVIII AT13148 perifosine",
 "PI3Ki": "pictilisib AZD6482 idelalisib TGX221 HS173 PIK294 PI3Ka_inh2 wortmannin LY294002 PI3Ki",
 "CK2i": "CX4945 TBCA",
 "MEKi": "selumetinib PD0325901 trametinib cobimetinib refametinib U0126 MEKi PD184352",
 "ERKi": "ERKi_SCH772984 ERKi_generic ravoxertinib ulixertinib",
 "RAFi": "vemurafenib dabrafenib SB590885 BRAFi LY3009120",
 "EGFRi": "gefitinib erlotinib lapatinib afatinib AG1478 AZD3759 PD153035 PD168393 CUDC101",
 "CDKi": "RO3306 dinaciclib AZD5438 flavopiridol",
}
A2C = {a: c for c, s in CLS.items() for a in s.split()}
