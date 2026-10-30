"""Shared helpers for I26 (perturbation audit of PSP kinase-substrate edges)."""
import numpy as np, pandas as pd

C26 = "/workspace/cache/b_I26"
R26 = "/workspace/results/b_I26"


def load_base():
    E = pd.read_parquet(f"{C26}/psp_edges_annot.parquet")
    Z = pd.read_parquet(f"{C26}/Z_lof.parquet").set_index("site")
    GC = pd.read_csv(f"{C26}/lof_conditions.tsv", sep="\t", dtype={"study": str})
    return E, Z, GC


def gold_sets(E, min_lt=3):
    """Per group: PSP in-vivo edges with site LT_LIT>=min_lt (fallback: all in-vivo edges if <5)."""
    out = {}
    for g, e in E[E.group.notna()].groupby("group"):
        s = e[e.in_vivo & (e.LT_LIT >= min_lt)].site.unique()
        if len(s) < 5:
            s = e[e.in_vivo].site.unique()
        out[g] = set(s)
    return out


def worked_table(Z, GC, gold, thr=-2.5, min_n=3):
    """Per (group, condition): control Stouffer z over gold sites; worked flag."""
    rows = []
    for g, d in GC.groupby("group"):
        cs = list(d.conditionID.unique())
        zz = Z.reindex(sorted(gold.get(g, ())))[cs]
        n = zz.notna().sum(); s = zz.sum()
        for c in cs:
            st = s[c] / np.sqrt(n[c]) if n[c] > 0 else np.nan
            rows.append((g, c, int(n[c]), float(s[c]), st, bool(n[c] >= min_n and st <= thr) if thr is not None else True))
    W = pd.DataFrame(rows, columns=["group", "conditionID", "ctrl_n", "ctrl_sum", "ctrl_stouffer", "worked"])
    return W.merge(GC[["group", "conditionID", "study", "agent"]], on=["group", "conditionID"])


def study_long(Z, GC, gold, thr=-2.5, min_n=3, sites=None):
    """Long table group, site, study, z (mean z over worked conditions of that study), ncond.
    Leave-one-out control for gold sites (their own z is removed from the control score)."""
    W = worked_table(Z, GC, gold, thr, min_n)
    out = []
    Zs = Z if sites is None else Z.reindex(sorted(set(sites) & set(Z.index)))
    for g, w in W.groupby("group"):
        for st, ws in w.groupby("study"):
            wc = list(ws.conditionID[ws.worked])
            if wc:
                sub = Zs[wc]
                m = sub.mean(axis=1); n = sub.notna().sum(axis=1)
                df = pd.DataFrame({"z": m, "ncond": n})
            else:
                df = pd.DataFrame({"z": pd.Series(dtype=float), "ncond": pd.Series(dtype=int)})
            # LOO for gold sites
            gs = [x for x in gold.get(g, ()) if x in Zs.index]
            if gs and thr is not None:
                allc = list(ws.conditionID)
                zg = Zs.loc[gs, allc]
                cs_ = ws.set_index("conditionID")
                n0 = cs_.ctrl_n.reindex(allc).values; s0 = cs_.ctrl_sum.reindex(allc).values
                v = zg.values
                meas = ~np.isnan(v)
                n1 = n0[None, :] - meas; s1 = s0[None, :] - np.where(meas, v, 0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    st1 = s1 / np.sqrt(n1)
                ok = (n1 >= min_n) & (st1 <= thr) & meas
                with np.errstate(invalid="ignore"):
                    zm = np.where(ok, v, 0).sum(1) / ok.sum(1)
                df = df.drop(index=[x for x in gs if x in df.index])
                df = pd.concat([df, pd.DataFrame({"z": zm, "ncond": ok.sum(1)}, index=gs)])
            df = df[df.ncond > 0]
            df.index.name = "site"
            out.append(df.reset_index().assign(group=g, study=st))
    L = pd.concat(out, ignore_index=True)
    return L, W


def classify(L, dec=-1.0, mean_unsup=-0.5, frac_unsup=0.25, studies_excl=None, studies_only=None):
    """Classify (group, site) from per-study z. Returns DataFrame with n_meas, n_dec, mean_z, cls."""
    if studies_excl is not None:
        L = L[~L.study.isin(studies_excl)]
    if studies_only is not None:
        L = L[L.study.isin(studies_only)]
    g = L.assign(d=(L.z <= dec)).groupby(["group", "site"], observed=True)
    A = g.agg(n_meas=("z", "size"), n_dec=("d", "sum"), mean_z=("z", "mean"), min_z=("z", "min")).reset_index()
    f = A.n_dec / A.n_meas
    A["cls"] = np.select([A.n_meas < 2, (A.n_dec >= 2) & (f >= 0.5), (f < frac_unsup) & (A.mean_z > mean_unsup)],
                         ["untestable", "supported", "unsupported"], "ambiguous")
    return A
