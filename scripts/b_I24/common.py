"""Shared definitions for b_I24 (signal-transmission gain)."""
import warnings; warnings.filterwarnings("ignore")
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
OUT = '/workspace/results/b_I24'
COHORTS = cl.COHORTS
U = dict(EGFR='P00533', ERBB2='P04626', ERBB3='P21860', MET='P08581', SHC1='P29353', GAB1='Q13480', PTPN11='Q06124',
         MAPK3='P27361', MAPK1='P28482', RPS6KA1='Q15418', RPS6KA3='P51812', AKT1S1='Q96B36', GSK3B='P49841', GSK3A='P49840',
         FOXO3='O43524', FOXO1='Q12778', TSC2='P49815', EIF4EBP1='Q13541', RPS6='P62753', MAPK9='P45984', JUN='P05412',
         ATF2='P15336', MAPK14='Q16539', MAPKAPK2='P49137', HSPB1='P04792', AKT2='P31751')
def S(g, *sites): return [f'{U[g]}_{s}' for s in sites]
STAGES = {
  'RTK':    S('EGFR', 'Y1092', 'Y1172', 'Y1197', 'Y1110', 'Y1016') + S('ERBB2', 'Y1248', 'Y877') + S('MET', 'Y1234', 'Y1235'),
  'EGFRpY': S('EGFR', 'Y1092', 'Y1172', 'Y1197', 'Y1110', 'Y1016'),
  'ADAPT':  S('SHC1', 'Y349', 'Y350', 'Y427') + S('GAB1', 'Y627', 'Y659') + S('PTPN11', 'Y542', 'Y580'),
  'ERK':    S('MAPK3', 'T202', 'Y204') + S('MAPK1', 'T185', 'Y187'),
  'RSK':    S('RPS6KA1', 'T359', 'S363', 'S380') + S('RPS6KA3', 'T365', 'S369', 'S386'),
  'AKTSUB': S('AKT1S1', 'T246') + S('GSK3B', 'S9') + S('GSK3A', 'S21') + S('FOXO3', 'S253') + S('FOXO1', 'S256') + S('TSC2', 'S939', 'S981', 'T1462'),
  'MTORC1': S('EIF4EBP1', 'T37', 'T46', 'S65', 'T70'),
  'S6':     S('RPS6', 'S235', 'S236', 'S240', 'S244'),
  'JNK':    S('MAPK9', 'T183', 'Y185'),
  'JNKSUB': S('JUN', 'S63', 'S73') + S('ATF2', 'T69', 'T71'),
  'P38':    S('MAPK14', 'T180', 'Y182'),
  'P38SUB': S('MAPKAPK2', 'T334') + S('HSPB1', 'S78', 'S82'),
}
# ERKSUB (PSP ERK1/2 substrate footprint) built in 01_stages.py
EDGES = [('RTK', 'ADAPT'), ('ADAPT', 'ERK'), ('ERK', 'RSK'), ('ERK', 'ERKSUB'), ('RTK', 'ERK'),
         ('ADAPT', 'AKTSUB'), ('AKTSUB', 'MTORC1'), ('MTORC1', 'S6'), ('JNK', 'JNKSUB'), ('P38', 'P38SUB')]
def ename(e): return f'{e[0]}>{e[1]}'

def zs(x):
    x = np.asarray(x, float); return (x - np.nanmean(x)) / np.nanstd(x)

def load_stages(variant='tc_padj'):
    return pd.read_parquet(f'{OUT}/tables/stages_{variant}.parquet')
import statsmodels.api as sm

def rlm_fit(y, X):
    """Huber robust regression; returns params, bse, resid (pandas)."""
    r = sm.RLM(y, sm.add_constant(X, has_constant='add'), M=sm.robust.norms.HuberT()).fit()
    return r

def meta_re(b, se):
    """DerSimonian-Laird random-effects meta. returns est, se, z, p, I2, tau2"""
    from scipy import stats
    b, se = np.asarray(b, float), np.asarray(se, float)
    ok = np.isfinite(b) & np.isfinite(se) & (se > 0); b, se = b[ok], se[ok]
    k = len(b)
    if k == 0: return dict(est=np.nan, se=np.nan, z=np.nan, p=np.nan, I2=np.nan, tau2=np.nan, k=0)
    w = 1 / se**2; fe = np.sum(w * b) / w.sum(); Q = np.sum(w * (b - fe)**2)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - np.sum(w**2) / w.sum())) if k > 1 else 0
    wr = 1 / (se**2 + tau2); est = np.sum(wr * b) / wr.sum(); s = np.sqrt(1 / wr.sum())
    z = est / s; p = 2 * stats.norm.sf(abs(z)); I2 = max(0, (Q - (k - 1)) / Q) if Q > 0 else 0
    return dict(est=est, se=s, z=z, p=p, I2=I2, tau2=tau2, k=k)
