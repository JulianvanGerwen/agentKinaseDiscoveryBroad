import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/b_I22'); from common_i22 import get
CA = '/workspace/cache/b_I25/'
DEF_STUDIES = {'21712546', '21857030', '26055452', '26270265', '27134283', '29970458'}
CC_CLASSES = {'CDK', 'CDK4_6', 'AURK_PLK'}

def load_test(kind='inh', include_cc=False):
    A, E = get(kind, st_only=False)
    keep = ~E.study.isin(DEF_STUDIES)
    if kind == 'inh' and not include_cc: keep &= ~E.cls.isin(CC_CLASSES)
    E = E[keep].reset_index(drop=True); A = A[E.exp.tolist()]
    A = A[A.notna().any(axis=1)]
    return A, E

def axes(index):
    return pd.read_parquet(CA + 'axes.parquet').reindex(index)
