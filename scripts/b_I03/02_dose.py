"""I03 step 2: dose coherence per kinase x site (decryptM curve fits from b_I17 + Klaeger Kinobeads Kd_app).
Experiments: single-agent kinase-inhibitor dose series with Klaeger profile. For kinase K and site s:
x_e = pKd_K(drug_e) (censored at 4 = not bound), y_e = site pEC50 if significant DECREASE (q<=.05, EC50 in range) else 4.
dose_r = Pearson(x, y) over experiments where s was fitted (>=5 exps, >=2 exps with pKd_K>=6 and >=2 with pKd_K<6).
dose_hit = # exps with pKd_K>=6.5 where s is a significant decrease; dose_nstrong = # such exps where s measured.
dose_gap = median(pEC50 - pKd_K) over hits (direct substrates expected ~ -1 to -3: ATP competition).
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
OUT = "/workspace/cache/b_I03"
F = pd.read_parquet("/workspace/results/b_I17/curves_all.parquet",
                    columns=["exp_id", "site", "delta", "sig", "ec_in_range", "pEC50"])
E = pd.read_csv("/workspace/cache/b_I17/experiments.tsv", sep="\t")
DM = pd.read_csv("/workspace/results/b_I17/drug_map.tsv", sep="\t").set_index("agent")
KD = pd.read_parquet("/workspace/cache/b_I17/klaeger_kd.parquet")
KC = pd.read_csv(f"{OUT}/kinase_conditions.tsv", sep="\t")
panel = sorted(KC.kinase.unique())

E = E[(E.pclass == "kinase_inhibitor") & E.agent.notna() & ~E.agent.str.contains(";", na=True)]
E = E[E.agent.isin(DM.index[DM.all_mapped])]
E = E[E.exp_id.isin(set(F.exp_id))]
exps = list(E.exp_id); drug = dict(zip(E.exp_id, E.agent.map(DM.klaeger)))
print(len(exps), "experiments", sorted(set(drug.values())))

# per-gene min Kd in Klaeger (rows may be 'A;B')
kd = {}
for name, row in KD.iterrows():
    for g in str(name).split(";"):
        kd[g] = np.fmin(kd.get(g, np.full(KD.shape[1], np.inf)), row.values.astype(float))
KDg = pd.DataFrame(kd, index=KD.columns).T

F = F[F.exp_id.isin(exps)]
F["y"] = np.where(F.sig & F.ec_in_range & (F.delta < 0), F.pEC50, 4.0)
Y = F.pivot_table(index="site", columns="exp_id", values="y", aggfunc="max")[exps]
H = F.assign(h=(F.sig & F.ec_in_range & (F.delta < 0)).astype(float)).pivot_table(index="site", columns="exp_id", values="h", aggfunc="max").reindex(index=Y.index, columns=exps)
Yv = Y.values; M = ~np.isnan(Yv); Hv = np.nan_to_num(H.values)
P = F.assign(p=np.where(F.sig & F.ec_in_range & (F.delta < 0), F.pEC50, np.nan)).pivot_table(index="site", columns="exp_id", values="p", aggfunc="max").reindex(index=Y.index, columns=exps).values
print("Y", Y.shape)
out = []
for k in panel:
    if k not in KDg.index:
        continue
    x = np.array([9 - np.log10(KDg.loc[k, drug[e]]) if np.isfinite(KDg.loc[k, drug[e]]) else 4.0 for e in exps])
    x = np.maximum(x, 4.0)
    X = np.where(M, x[None, :], np.nan)
    n = M.sum(1)
    hi = (M & (x[None, :] >= 6)).sum(1); lo = (M & (x[None, :] < 6)).sum(1)
    mx = np.nanmean(X, 1, keepdims=True); my = np.nanmean(Yv, 1, keepdims=True)
    cov = np.nansum((X - mx) * (Yv - my), 1)
    vx = np.nansum((X - mx) ** 2, 1); vy = np.nansum((Yv - my) ** 2, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = cov / np.sqrt(vx * vy)
    ok = (n >= 5) & (hi >= 2) & (lo >= 2)
    r = np.where(ok & (vy > 0), r, np.where(ok, 0.0, np.nan))
    strong = x >= 6.5
    nstrong = (M & strong[None, :]).sum(1); nhit = (Hv[:, strong] > 0).sum(1)
    gap = np.nanmedian(np.where(Hv[:, strong] > 0, P[:, strong] - x[None, strong], np.nan), 1) if strong.any() else np.full(len(Y), np.nan)
    df = pd.DataFrame(dict(kinase=k, site=Y.index.values, dose_r=r.astype(np.float32), dose_n=n, dose_nstrong=nstrong,
                           dose_hit=nhit, dose_gap=gap))
    df = df[ok | (nstrong > 0)]
    print(k, "strong drugs:", sorted({drug[e] for e, s in zip(exps, strong) if s}), "evaluable sites:", int(ok.sum()))
    out.append(df)
D = pd.concat(out, ignore_index=True)
D.to_parquet(f"{OUT}/dose_long.parquet"); print(D.shape)
