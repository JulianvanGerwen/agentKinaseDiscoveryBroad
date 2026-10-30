"""Signed-set pair scoring with exclusive (non-shared, non-same-protein) site sets, vectorised over samples/experiments."""
import numpy as np, pandas as pd

def prepare_sets(sets, universe):
    """sets: dict name -> Series(sign, index=site). Restrict to universe."""
    out = {}
    for k, s in sets.items():
        s = s[~s.index.duplicated()]; s = s[s.index.isin(universe)]
        out[k] = s
    return out

def pair_scores(Z, sets, pairs, minsites=5):
    """Z: site x sample DataFrame (standardised; NaN allowed). For each pair (a,b): score_a on sites of a not in b and not on a protein
    with sites in b, and vice versa. Returns dict (a,b) -> (sa, sb) arrays over samples (NaN when < minsites)."""
    idx = {s: i for i, s in enumerate(Z.index)}; X0 = np.nan_to_num(Z.values); M = np.isfinite(Z.values).astype(float)
    prot = {k: set(x.split('_')[0] for x in s.index) for k, s in sets.items()}
    out = {}
    for a, b in pairs:
        sa, sb = sets[a], sets[b]
        ea = sa[[x.split('_')[0] not in prot[b] for x in sa.index]]
        eb = sb[[x.split('_')[0] not in prot[a] for x in sb.index]]
        res = []
        for e in (ea, eb):
            if len(e) < minsites: res.append(None); continue
            ii = np.array([idx[x] for x in e.index]); w = e.values[:, None]
            num = (X0[ii] * w).sum(0); n = M[ii].sum(0)
            with np.errstate(all='ignore'): v = num / n
            v[n < minsites] = np.nan; res.append(v)
        out[(a, b)] = tuple(res)
    return out

def rcorr(x, y, minn=10):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < minn: return np.nan, int(m.sum())
    return np.corrcoef(x[m], y[m])[0, 1], int(m.sum())
