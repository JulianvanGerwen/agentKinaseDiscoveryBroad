"""Characterise sites driving PPP2R2A-protein association (pc5 model): per-site meta across cohorts; GO BP enrichment of top up-sites' proteins."""
import gzip, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
R = "/workspace/results/b_I13"; TAG = "_pc5"
F = pd.read_parquet(f"{R}/cache/site_features.parquet").set_index("site")
B = {c: pd.read_parquet(f"{R}/cache/stage1_{c}{TAG}.parquet") for c in cl.COHORTS}
for ex in ["PPP2R2A_protlow", "PPP2R2A_cnvloss"]:
    M = pd.DataFrame({c: b[ex] for c, b in B.items() if ex in b})
    k = M.notna().sum(1); M = M[k >= 5]
    # per-cohort z of site beta (beta / sd of betas across sites in cohort)
    Z = M / M.std()
    meta = pd.DataFrame({"k": M.notna().sum(1), "mean_beta": M.mean(1), "stouffer": Z.sum(1) / np.sqrt(Z.notna().sum(1)), "npos": (M > 0).sum(1)})
    meta = meta.join(F[["residue", "window", "uniprot", "B55_lit", "bipartite_basic", "isT", "P1"]])
    meta["p"] = 2 * stats.norm.sf(meta.stouffer.abs())
    from statsmodels.stats.multitest import multipletests
    meta["fdr"] = multipletests(meta.p, method="fdr_bh")[1]
    meta.sort_values("stouffer", ascending=False).to_csv(f"{R}/tables/t3_sites_meta_{ex}.tsv", sep="\t")
    print(ex, len(meta), "FDR<0.1 up/down:", ((meta.fdr < 0.1) & (meta.stouffer > 0)).sum(), ((meta.fdr < 0.1) & (meta.stouffer < 0)).sum())
    st = meta[meta.residue.isin(["S", "T"])]
    for f in ["bipartite_basic", "isT", "P1"]:
        top = st.stouffer > st.stouffer.quantile(0.95); bot = st.stouffer < st.stouffer.quantile(0.05)
        print(f"  {f}: top5% {st[top][f].mean():.3f}  bottom5% {st[bot][f].mean():.3f}  all {st[f].mean():.3f}")
# GO BP enrichment for proteins with a top-5% up site (PPP2R2A_protlow) vs all measured proteins
meta = pd.read_csv(f"{R}/tables/t3_sites_meta_PPP2R2A_protlow.tsv", sep="\t", index_col=0)
obo = {}; cur = None
for line in open("<HACKATHON_ROOT>/kinaseSubMain/external/go/go-basic.obo"):
    if line.startswith("id: GO:"): cur = line.split()[1]
    elif line.startswith("name: ") and cur: obo[cur] = line[6:].strip()
    elif line.startswith("namespace: ") and cur: obo[cur] = (obo[cur], line.split()[1])
g2p = {}
for line in gzip.open("<HACKATHON_ROOT>/kinaseSubMain/external/go/goa_human.gaf.gz", "rt"):
    if line.startswith("!"): continue
    f = line.split("\t")
    if f[8] != "P" or "NOT" in f[3]: continue
    g2p.setdefault(f[4], set()).add(f[1])
univ = set(meta.uniprot.dropna()); hit = set(meta[meta.stouffer > meta.stouffer.quantile(0.95)].uniprot.dropna())
low = set(meta[meta.stouffer < meta.stouffer.quantile(0.05)].uniprot.dropna())
rows = []
for go, ps in g2p.items():
    ps = ps & univ
    if not 15 <= len(ps) <= 500: continue
    for lab, H in [("up", hit), ("down", low)]:
        a = len(ps & H); b = len(H) - a; c_ = len(ps) - a; d = len(univ) - a - b - c_
        orr, p = stats.fisher_exact([[a, b], [c_, d]], alternative="greater")
        rows.append(dict(direction=lab, go=go, term=obo.get(go, ("?",))[0], n_term=len(ps), n_hit=a, odds=orr, p=p))
G = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
G["fdr"] = G.groupby("direction").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
G.sort_values("p").to_csv(f"{R}/tables/t3_GO_PPP2R2A_protlow_sites.tsv", sep="\t", index=False)
for lab in ["up", "down"]: print(lab); print(G[G.direction == lab].sort_values("p").head(10)[["term", "n_term", "n_hit", "odds", "p", "fdr"]].round(4).to_string(index=False))
