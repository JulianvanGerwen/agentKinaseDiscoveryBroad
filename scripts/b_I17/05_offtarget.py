"""Functional off-targets.
(A) Within-drug potency matching: for each experiment (single Klaeger-profiled drug) and each Klaeger kinase k
    (Kd <= 10 uM, not in the on-target cluster T), sites whose pEC50 is within 0.5 log of pKd_k are tested for
    enrichment of k among their kinaseSub top5 predicted kinases (Fisher, vs all other regulated down-sites).
    Control: drug-label shuffle = use kinase profile of every other drug on this experiment's sites.
(B) Cross-drug regression (A549 panel + all lines pooled): per site, pEC50 (censored to 4.5 if not regulated)
    vs pKd of each kinase across drugs; leave-one-drug-out (LODO) prediction; null = drug-label shuffle.
Outputs: results/b_I17/offtarget_matching.tsv, offtarget_matching_null.tsv, offtarget_crossdrug.tsv
"""
import os
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, spearmanr
OUT = "/workspace/cache/b_I17"; RES = "/workspace/results/b_I17"
F = pd.read_parquet(f"{RES}/curves_all.parquet")
E = pd.read_csv(f"{OUT}/experiments.tsv", sep="\t").set_index("exp_id")
K = pd.read_parquet(f"{OUT}/klaeger_kd.parquet")
DM = pd.read_csv(f"{RES}/drug_map.tsv", sep="\t").set_index("agent")
A = pd.read_csv(f"{RES}/shift_sites.tsv", sep="\t")
P5 = pd.read_parquet(f"{OUT}/pred_top5.parquet")
pred = P5[P5["rank"] <= 5].groupby("site_id").gene.apply(set)
pred_univ = set(P5.gene)
rng = np.random.default_rng(1)

exps = [ex for ex in E.index if E.loc[ex, "pclass"] == "kinase_inhibitor" and isinstance(E.loc[ex, "agent"], str)
        and ";" not in E.loc[ex, "agent"] and E.loc[ex, "agent"] in DM.index and DM.loc[E.loc[ex, "agent"], "all_mapped"]
        and ex in set(F.exp_id)]
def profile(kd):
    c = K[kd]; c = c[np.isfinite(c)].dropna(); c = c[c <= 10000]
    out = {}
    for name, v in c.items():
        for gname in str(name).split(";"):
            if gname in pred_univ: out[gname] = min(v, out.get(gname, np.inf))
    return pd.Series(out, dtype=float)
profiles = {DM.loc[E.loc[ex, "agent"], "klaeger"]: profile(DM.loc[E.loc[ex, "agent"], "klaeger"]) for ex in exps}


def match_test(sites, prof, exclude):
    res = []
    if len(sites) == 0: return res
    has = np.array([s in pred.index for s in sites.site], dtype=bool)
    sites = sites.loc[has]
    pk = np.array([pred[s] for s in sites.site], dtype=object)
    for k, kd in prof.items():
        if k in exclude: continue
        pkd = 9 - np.log10(kd)
        m = (sites.pEC50 - pkd).abs().values <= 0.5
        hit = np.array([k in s for s in pk])
        a, b, c, d = (m & hit).sum(), (m & ~hit).sum(), (~m & hit).sum(), (~m & ~hit).sum()
        if m.sum() < 5 or hit.sum() < 3: continue
        orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        res.append(dict(kinase=k, kd_nM=kd, pKd=pkd, n_match=int(m.sum()), n_hit_match=int(a), n_hit_other=int(c),
                        frac_hit_match=a / max(m.sum(), 1), frac_hit_other=c / max((~m).sum(), 1),
                        OR=(a + .5) * (d + .5) / ((b + .5) * (c + .5)), p=p))
    return res


rows, nrows = [], []
for ex in exps:
    agent = E.loc[ex, "agent"]; kd = DM.loc[agent, "klaeger"]
    sig = F[(F.exp_id == ex) & F.sig & (F.delta < 0)]
    T = set(A.loc[A.exp_id == ex, "target"]) if (A.exp_id == ex).any() else set()
    # exclude on-target cluster: kinases within 3x of the best designated target (from shift table) + designated
    prof = profiles[kd]
    if (A.exp_id == ex).any():
        pkd_t = A.loc[A.exp_id == ex, "pKd_target"].iloc[0]
        T |= set(prof[(9 - np.log10(prof)) >= pkd_t - np.log10(3)].index)
    for r in match_test(sig, prof, T):
        rows.append(dict(exp_id=ex, agent=agent, cell_line=E.loc[ex, "cell_line"], **r))
    for kd2, prof2 in profiles.items():
        if kd2 == kd: continue
        for r in match_test(sig, prof2, set()):
            nrows.append(dict(exp_id=ex, agent=agent, other_drug=kd2, **r))
R = pd.DataFrame(rows); N = pd.DataFrame(nrows)
from statsmodels.stats.multitest import multipletests
if len(R):
    R["q"] = multipletests(R.p, method="fdr_bh")[1]
R.sort_values("p").to_csv(f"{RES}/offtarget_matching.tsv", sep="\t", index=False)
N.to_csv(f"{RES}/offtarget_matching_null.tsv", sep="\t", index=False)
print(R.sort_values("p").head(30).round(4).to_string())
print("true-profile tests:", len(R), "frac p<0.01:", (R.p < 0.01).mean() if len(R) else None,
      "| shuffled-profile tests:", len(N), "frac p<0.01:", (N.p < 0.01).mean() if len(N) else None)

# (B) cross-drug regression
rowsB = []
for scope in ["A549", "all"]:
    ex_s = [ex for ex in exps if scope == "all" or E.loc[ex, "cell_line"] == scope]
    # one value per drug (average replicate/cell-line experiments)
    Y = {}
    for ex in ex_s:
        d = DM.loc[E.loc[ex, "agent"], "klaeger"]
        f = F[F.exp_id == ex].set_index("site")
        y = pd.Series(4.5, index=f.index)
        dn = f.sig & (f.delta < 0)
        y[dn] = f.pEC50[dn]
        Y.setdefault(d, []).append(y)
    Y = pd.DataFrame({d: pd.concat(v, axis=1).mean(axis=1) for d, v in Y.items()})
    drugs = list(Y.columns)
    Y = Y[(Y.notna().sum(axis=1) >= 6) & ((Y > 4.5).sum(axis=1) >= 2)]
    KP = pd.DataFrame({d: K[d] for d in drugs})
    KP.index = [str(i).split(";")[0] for i in KP.index]
    KP = KP[~KP.index.duplicated()]
    KP = (9 - np.log10(KP.replace(np.inf, 10 ** 4.5))).where(KP.notna())
    KP = KP[(KP > 5).sum(axis=1) >= 2]
    KP = KP[KP.notna().sum(axis=1) >= 6]
    print(scope, "drugs", len(drugs), "sites", len(Y), "kinases", len(KP))
    Yv = Y.values; Kv = KP.values

    def best_corr(Yv, Kv):
        out_k = np.full(len(Yv), -1); out_r = np.full(len(Yv), np.nan)
        for i in range(len(Yv)):
            m0 = np.isfinite(Yv[i])
            best = (-2, -1)
            for j in range(len(Kv)):
                m = m0 & np.isfinite(Kv[j])
                if m.sum() < 6 or np.std(Kv[j, m]) == 0 or np.std(Yv[i, m]) == 0: continue
                r = np.corrcoef(Yv[i, m], Kv[j, m])[0, 1]
                if r > best[0]: best = (r, j)
            out_r[i], out_k[i] = best
        return out_k, out_r
    bk, br = best_corr(Yv, Kv)
    # LODO: for each site, for each held-out drug, choose best kinase on remaining drugs, predict by linear fit
    def lodo(Yv, Kv):
        err = []
        for i in range(len(Yv)):
            for h in range(Yv.shape[1]):
                if not np.isfinite(Yv[i, h]): continue
                tr = np.isfinite(Yv[i]); tr[h] = False
                best = (-2, -1)
                for j in range(len(Kv)):
                    m = tr & np.isfinite(Kv[j])
                    if m.sum() < 5 or not np.isfinite(Kv[j, h]) or np.std(Kv[j, m]) == 0 or np.std(Yv[i, m]) == 0: continue
                    r = np.corrcoef(Yv[i, m], Kv[j, m])[0, 1]
                    if r > best[0]: best = (r, j)
                if best[1] < 0: continue
                j = best[1]; m = tr & np.isfinite(Kv[j])
                b = np.polyfit(Kv[j, m], Yv[i, m], 1)
                err.append((i, h, abs(np.polyval(b, Kv[j, h]) - Yv[i, h]), abs(Yv[i, tr].mean() - Yv[i, h])))
        return pd.DataFrame(err, columns=["i", "h", "err_kinase", "err_mean"])
    sub = rng.choice(len(Yv), size=min(1500, len(Yv)), replace=False)
    L = lodo(Yv[sub], Kv)
    nullL = []
    for p in range(5):
        perm = rng.permutation(Kv.shape[1])
        Ln = lodo(Yv[sub], Kv[:, perm]); nullL.append(Ln.err_kinase.median())
    rowsB.append(dict(scope=scope, n_drugs=len(drugs), n_sites=len(Y), n_kinases=len(KP),
                      lodo_median_err_kinase=L.err_kinase.median(), lodo_median_err_mean=L.err_mean.median(),
                      frac_kinase_better=(L.err_kinase < L.err_mean).mean(), null_shuffle_median_err=np.median(nullL),
                      null_shuffle_runs=";".join(f"{v:.3f}" for v in nullL)))
    T = pd.DataFrame(dict(site=Y.index, best_kinase=[KP.index[j] if j >= 0 else None for j in bk], r=br))
    T["pred_top5"] = [";".join(sorted(pred.get(s, set()))) for s in T.site]
    T["motif_agrees"] = [bk_ in pred.get(s, set()) if bk_ else False for s, bk_ in zip(T.site, T.best_kinase)]
    T["scope"] = scope
    T.to_csv(f"{RES}/offtarget_crossdrug_{scope}.tsv", sep="\t", index=False)
B = pd.DataFrame(rowsB); B.to_csv(f"{RES}/offtarget_crossdrug.tsv", sep="\t", index=False)
print(B.to_string())
