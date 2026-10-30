"""P1 (pre-registered) + S1/S2/S5 on Hijazi. Out: tables/T3_hijazi_pairs.tsv, T4_primary_results.tsv"""
import numpy as np, pandas as pd
T = "/workspace/results/b_I23/tables/"
F = pd.read_csv(T + "T1_features.tsv", sep="\t"); S = pd.read_csv(T + "T2_sens_matched.tsv", sep="\t")
H = F[F.study == 31959955].copy() if F.study.dtype != object else F[F.study == "31959955"].copy()
Sh = S[S.cell.isin(["HL60", "MCF7", "NTERA2"])]
pick = {}
for d, g in Sh.groupby("drug"):
    n2 = (g.source == "GDSC2").sum(); n1 = (g.source == "GDSC1").sum()
    pick[d] = "GDSC2" if n2 >= 2 else ("GDSC1" if n1 >= 2 else None)
Sh = Sh[Sh.apply(lambda r: pick.get(r.drug) == r.source, axis=1)]
D = H.merge(Sh, on=["drug", "cell"])
D = D[D.groupby("drug").cell.transform("size") >= 2]
D["logfu"] = np.log(D.f_up + .002); D["logfd"] = np.log(D.f_down + .002)
D.to_csv(T + "T3_hijazi_pairs.tsv", sep="\t", index=False)
print("drugs", D.drug.nunique(), "pairs", len(D), D.groupby("drug").size().value_counts().to_dict(), D.source.value_counts().to_dict())

_MC = {}
def _M(df):
    key = tuple(df.drug) + tuple(df.cell)
    if key not in _MC:
        Dm = pd.get_dummies(df.drug).values.astype(float); Cm = pd.get_dummies(df.cell).values.astype(float)[:, 1:]
        A = np.column_stack([Dm, Cm]); _MC.clear(); _MC[key] = np.eye(len(df)) - A @ np.linalg.pinv(A)
    return _MC[key]
def resid(df, col):  # exact two-way FE residualisation
    return _M(df) @ df[col].astype(float).values

def fe_beta(df, y, xs):
    Y = resid(df, y); Xm = np.column_stack([resid(df, x) for x in xs])
    b, *_ = np.linalg.lstsq(Xm, Y, rcond=None)
    return b, Y, Xm

def run(df, y, xs, nperm=10000, nboot=2000, seed=0, permcol=None):
    df = df.dropna(subset=[y] + xs).reset_index(drop=True)
    b, Y, Xm = fe_beta(df, y, xs)
    rng = np.random.default_rng(seed)
    # standardized beta for first x: per SD of residualised x, outcome per SD of residualised y
    sb = b[0] * Xm[:, 0].std() / Y.std()
    # permutation within drug of first x
    xcol = permcol or xs[0]; idx = df.groupby("drug").indices
    null = []
    for _ in range(nperm):
        perm = np.arange(len(df))
        for d, ii in idx.items(): perm[ii] = rng.permutation(ii)
        d2 = df.copy(); d2[xcol] = df[xcol].values[perm]
        null.append(fe_beta(d2, y, xs)[0][0])
    null = np.array(null); p = (np.sum(np.abs(null) >= abs(b[0])) + 1) / (nperm + 1)
    drugs = df.drug.unique(); bs = []
    for _ in range(nboot):
        dd = rng.choice(drugs, len(drugs)); parts = []
        for i, d in enumerate(dd):
            q = df[df.drug == d].copy(); q["drug"] = f"{d}_{i}"; parts.append(q)
        q = pd.concat(parts); 
        try: bb, YY, XX = fe_beta(q, y, xs); bs.append(bb[0] * XX[:, 0].std() / YY.std())
        except Exception: pass
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(y=y, x=xs[0], covars=";".join(xs[1:]), n=len(df), n_drugs=df.drug.nunique(), beta=b[0], std_beta=sb,
                ci_lo=lo, ci_hi=hi, p_perm=p)

if __name__ == "__main__":
    R = []
    R.append(dict(test="P1_primary", **run(D, "auc", ["LR"])))
    R.append(dict(test="P1_lnIC50", **run(D, "ln_ic50", ["LR"], nperm=2000, nboot=500)))
    R.append(dict(test="S1_logfu_adj_fd", **run(D, "auc", ["logfu", "logfd"], nperm=2000, nboot=500)))
    R.append(dict(test="S1_logfd_adj_fu", **run(D, "auc", ["logfd", "logfu"], nperm=2000, nboot=500)))
    R.append(dict(test="S5_primary_only_logfd", **run(D, "auc", ["logfd"], nperm=2000, nboot=500)))
    R.append(dict(test="S1_thr0.58", **run(D.assign(LR58=np.log((D.f_up58+.002)/(D.f_down58+.002))), "auc", ["LR58"], nperm=2000, nboot=500)))
    R.append(dict(test="S1_sumLR", **run(D.assign(sLR=np.log(D.sum_up/D.sum_down)), "auc", ["sLR"], nperm=2000, nboot=500)))
    R.append(dict(test="S2_LR_adj_targetz", **run(D, "auc", ["LR", "target_z"], nperm=2000, nboot=500)))
    R.append(dict(test="S2_targetz_alone", **run(D, "auc", ["target_z"], nperm=2000, nboot=500)))
    R.append(dict(test="S5_LR_halfsites", **run(D, "auc", ["LR_half"], nperm=2000, nboot=500)))
    R = pd.DataFrame(R); R.to_csv(T + "T4_primary_results.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); print(R.round(4).to_string())
