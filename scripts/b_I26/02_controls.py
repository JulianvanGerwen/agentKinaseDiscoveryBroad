"""I26 step 2: random-relabel control (same group, random measured sites matched on n_meas and residue class,
excluding the group's PSP sites) + threshold sensitivity grid.
Outputs tables/random_control.tsv, tables/sensitivity_grid.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *

E, Z, GC = load_base()
gold = gold_sets(E)
rng = np.random.default_rng(26)
PSPg = E[E.group.notna()][["group", "site"]].drop_duplicates()
psp_key = set(PSPg.group + "|" + PSPg.site)
res_y = lambda s: s.str.split("_").str[1].str[0].eq("Y")


def random_control(L, A_psp, nrep=20, **kw):
    """A_psp: classification of PSP edges (testable rows). Returns class fractions PSP vs random."""
    A_all = classify(L, **kw)
    A_all = A_all[A_all.n_meas >= 2]
    A_all = A_all[~(A_all.group.astype(str) + "|" + A_all.site).isin(psp_key)]
    A_all["Y"] = res_y(A_all.site)
    rows = []
    for g, ap in A_psp.groupby("group"):
        pool = A_all[A_all.group == g]
        ap = ap.assign(Y=res_y(ap.site), nm=ap.n_meas.clip(upper=6))
        pool = pool.assign(nm=pool.n_meas.clip(upper=6))
        pg = {k: v for k, v in pool.groupby(["Y", "nm"])}
        for (y, nm), a in ap.groupby(["Y", "nm"]):
            p = pg.get((y, nm))
            if p is None or len(p) < 5:
                continue
            idx = rng.integers(0, len(p), size=len(a) * nrep)
            r = p.cls.values[idx]
            for c in ["supported", "unsupported", "ambiguous"]:
                rows.append((g, y, nm, c, (a.cls == c).sum(), (r == c).sum() / nrep, len(a)))
    return pd.DataFrame(rows, columns=["group", "Y", "nm", "cls", "n_psp", "n_rand_exp", "n_edges"])


if __name__ == "__main__":
    grid = []
    out_rc = None
    for min_lt in [3, 1]:
        gold = gold_sets(E, min_lt)
        for thr in [None, -1.5, -2.5, -4.0]:
            L, W = study_long(Z, GC, gold, thr=thr)
            for dec in [-0.5, -1.0, -1.5]:
                A = classify(L, dec=dec)
                Ap = PSPg.merge(A, on=["group", "site"])
                Ap = Ap[Ap.n_meas >= 2]
                rc = random_control(L, Ap, dec=dec)
                tot = rc.groupby("cls")[["n_psp", "n_rand_exp"]].sum()
                n = rc.drop_duplicates(["group", "Y", "nm"]).n_edges.sum()
                r = dict(gold_min_LT=min_lt, ctrl_thr=thr if thr is not None else "none", dec_thr=dec, n_testable=n,
                         n_groups=Ap.group.nunique())
                for c in ["supported", "unsupported", "ambiguous"]:
                    r[f"psp_{c}"] = tot.loc[c, "n_psp"] / n; r[f"rand_{c}"] = tot.loc[c, "n_rand_exp"] / n
                grid.append(r); print(r)
                if min_lt == 3 and thr == -2.5 and dec == -1.0:
                    out_rc = rc
    G = pd.DataFrame(grid)
    G.to_csv(f"{R26}/tables/sensitivity_grid.tsv", sep="\t", index=False, float_format="%.3f")
    rcg = out_rc.groupby(["group", "cls"])[["n_psp", "n_rand_exp"]].sum().unstack()
    out_rc.to_csv(f"{R26}/tables/random_control_strata.tsv", sep="\t", index=False, float_format="%.3f")
    rcg.columns = [f"{a}_{b}" for a, b in rcg.columns]
    nn = out_rc.drop_duplicates(["group", "Y", "nm"]).groupby("group").n_edges.sum()
    rcg["n_edges_matched"] = nn
    rcg.to_csv(f"{R26}/tables/random_control.tsv", sep="\t", float_format="%.2f")
    print(rcg.round(1).to_string())
