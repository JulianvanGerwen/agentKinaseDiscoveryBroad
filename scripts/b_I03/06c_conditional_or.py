"""I03 step 6c: conditional value of new evidence among high-prior pairs (kinaseSub-DEV pct>=95).
Logit: y ~ evidence + C(ks_tier) + C(kinase) + log(n_st) [+ log site n_studies]; PSP labels (all proteins & TEST proteins) and
external curated non-PSP labels. Control: same model with evidence shuffled within kinase x n_st bin (coverage-matched random sites).
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd, statsmodels.formula.api as smf
OUT = "/workspace/cache/b_I03"; RES = "/workspace/results/b_I03/tables"
X = pd.read_parquet(f"{OUT}/scored.parquet")
X = X[X.ks_dev_pct >= 95].copy()
X["ks_tier"] = pd.cut(X.ks_dev_pct, [94, 99, 99.7, 100.1], labels=["95-99", "99-99.7", "gt99.7"]).astype(str)
X["drop2"] = (X.n_st_drop >= 2).astype(int); X["drop1"] = (X.n_st_drop == 1).astype(int)
X["tum_pos"] = (X.tum_z > 3).astype(int); X["tum_meas"] = X.tum_z.notna().astype(int)
X["dose_pos"] = (X.dose_r > 0.3).astype(int); X["dose_meas"] = X.dose_r.notna().astype(int)
X["lnst"] = np.log(X.n_st); X["lstud"] = np.log1p(X.n_studies.astype(float))
X["motif_hi"] = (X.motif >= 90).astype(int)
rng = np.random.default_rng(0)
rows = []
def run(D, y, label, shuffle=False):
    D = D.copy()
    if shuffle:
        D["nb"] = pd.qcut(D.n_st.rank(method="first"), 5, labels=False)
        for c in ["drop2", "drop1", "tum_pos", "tum_meas", "dose_pos", "dose_meas"]:
            D[c] = D.groupby(["kinase", "nb"])[c].transform(lambda v: v.sample(frac=1, random_state=int(rng.integers(1e9))).values)
    f = f"{y} ~ drop2 + drop1 + tum_pos + tum_meas + dose_pos + dose_meas + motif_hi + C(ks_tier) + lnst + lstud + C(kinase)"
    keep = D.groupby("kinase")[y].transform("sum") > 0
    m = smf.logit(f, D[keep]).fit(disp=0, method="bfgs", maxiter=2000)
    for t in ["drop2", "drop1", "tum_pos", "dose_pos", "motif_hi"]:
        rows.append(dict(label=label + (" [shuffled control]" if shuffle else ""), term=t, OR=np.exp(m.params[t]),
                         lo=np.exp(m.conf_int().loc[t, 0]), hi=np.exp(m.conf_int().loc[t, 1]), p=m.pvalues[t],
                         n=int(keep.sum()), n_pos=int(D[keep][y].sum())))
run(X, "y_psp", "PSP, all proteins")
run(X[X.test_protein == 1], "y_psp", "PSP, TEST proteins")
run(X[X.y_psp == 0], "y_ext", "external curated non-PSP")
for _ in range(1):
    run(X, "y_psp", "PSP, all proteins", shuffle=True)
    run(X[X.y_psp == 0], "y_ext", "external curated non-PSP", shuffle=True)
R = pd.DataFrame(rows); R.to_csv(f"{RES}/val_conditional_OR.tsv", sep="\t", index=False); print(R.round(4).to_string())
