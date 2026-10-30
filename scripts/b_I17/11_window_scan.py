"""Potency-resolved control for functional off-targets: enrichment (OR) of kinase-k predicted substrates among
down-sites in sliding pEC50 windows (width 1 log) - should peak near pKd_k (Kinobeads) if engagement is causal.
Output: results/b_I17/offtarget_window_scan.tsv, fig5_offtarget_window.png"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES = "/workspace/results/b_I17"; OUT = "/workspace/cache/b_I17"
F = pd.read_parquet(f"{RES}/curves_all.parquet")
P5 = pd.read_parquet(f"{OUT}/pred_top5.parquet"); pred = P5[P5["rank"] <= 5].groupby("site_id").gene.apply(set)
cases = [("DS0005_Afatinib0nM_0", "MAPKAPK2", 817), ("DS0006_Dasatinib0nM_0", "MAPK14", 945), ("DS0040_Dasatinib0nM_0", "MAPKAPK2", 538),
         ("DS0015_Dasatinib0nM_0", "MAPKAPK2", 538), ("DS0015_Dasatinib0nM_1", "MAPK14", 945), ("DS0006_Dasatinib0nM_0", "MAPKAPK2", 538),
         ("DS0005_Afatinib0nM_0", "EGFR", 2)]
rows = []
for ex, k, kd in cases:
    d = F[(F.exp_id == ex) & F.sig & (F.delta < 0)]
    d = d[d.site.isin(pred.index)]
    hit = np.array([k in pred[s] for s in d.site])
    for c in np.arange(4.5, 9.51, 0.25):
        m = (d.pEC50 - c).abs().values <= 0.5
        a, b, cc, dd = (m & hit).sum(), (m & ~hit).sum(), (~m & hit).sum(), (~m & ~hit).sum()
        rows.append(dict(exp_id=ex, kinase=k, pKd=9 - np.log10(kd), center=c, n=int(m.sum()), frac_hit=a / max(m.sum(), 1),
                         log2OR=np.log2((a + .5) * (dd + .5) / ((b + .5) * (cc + .5)))))
W = pd.DataFrame(rows); W.to_csv(f"{RES}/offtarget_window_scan.tsv", sep="\t", index=False)
fig, axs = plt.subplots(1, len(cases), figsize=(2.3 * len(cases), 2.6), sharey=True)
for ax, (ex, k, kd) in zip(axs, cases):
    w = W[(W.exp_id == ex) & (W.kinase == k) & (W.n >= 5)]
    ax.plot(w.center, w.log2OR, color="#2a78d6", lw=2); ax.axvline(9 - np.log10(kd), color="#eb6834", ls=":", lw=1.2)
    ax.axhline(0, color="#8a8984", lw=.6); ax.set_title(f"{ex.split('_')[1].replace('0nM','')} {ex.split('_')[0]}\n{k}", fontsize=7, loc="left")
    ax.set_xlabel("window centre pEC50", fontsize=7); ax.tick_params(labelsize=7)
    pk = w.loc[w.log2OR.idxmax()] if len(w) else None
    if pk is not None: print(ex, k, "pKd", round(9 - np.log10(kd), 2), "peak at", pk.center, "log2OR", round(pk.log2OR, 2))
axs[0].set_ylabel("log2 OR predicted substrate", fontsize=7)
fig.suptitle("Off-target substrate enrichment vs potency window (orange = Kinobeads pKd)", fontsize=8, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(f"{RES}/fig5_offtarget_window.png")
