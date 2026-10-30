import numpy as np, pandas as pd
from scipy import stats
CA = "/workspace/results/b_I19/cache/"
def experiments():
    E = pd.read_csv(CA + 'experiments.tsv', sep='\t', dtype={'study': str}).fillna({'binders': '', 'primary': '', 'klaeger1000': '', 'hijazi': ''})
    E['bset'] = E.binders.str.split(';').apply(lambda l: frozenset(x for x in l if x))
    E['grp'] = E.study + '|' + E.cell_line.astype(str)
    g = E.groupby('grp').agg(n=('exp', 'size'), ncls=('cls', 'nunique'))
    E['centerable'] = E.grp.map((g.n >= 5) & (g.ncls >= 3))
    return E
def center(M, E):
    """Within study x cell-line centering (subtract median over that group's inhibitor experiments) for groups with >=5 exps & >=3 classes;
    other experiments are left as is (vs-DMSO)."""
    M = M[E.exp].copy()
    for grp, d in E[E.centerable].groupby('grp'):
        cols = d.exp.tolist(); M[cols] = M[cols].sub(M[cols].median(axis=1), axis=0)
    return M
def dl_meta(y, v):
    """DerSimonian-Laird random effects. y, v: arrays (k,). returns est, se, p, tau2, I2"""
    k = len(y)
    w = 1 / v; fe = (w * y).sum() / w.sum()
    if k == 1: return fe, np.sqrt(1 / w.sum()), 2 * stats.norm.sf(abs(fe) / np.sqrt(1 / w.sum())), 0.0, 0.0
    Q = (w * (y - fe) ** 2).sum(); c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0, (Q - (k - 1)) / c); ws = 1 / (v + tau2)
    est = (ws * y).sum() / ws.sum(); se = np.sqrt(1 / ws.sum())
    I2 = max(0, (Q - (k - 1)) / Q) if Q > 0 else 0
    return est, se, 2 * stats.norm.sf(abs(est / se)), tau2, I2
