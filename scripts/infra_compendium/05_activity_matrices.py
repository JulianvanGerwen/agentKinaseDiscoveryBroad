"""Kinase x condition activity matrices.

Method A  core-v3 (cache/compendium/corev3_long.parquet from 03_score_corev3.py).
Method B  'psp_z': independent, simple. Library = PhosphoSitePlus Kinase_Substrate_Dataset only (kinase & substrate
          organism human/mouse/rat; non-human substrate sites mapped to the human site via SITE_GRP_ID; kinase mapped to
          the 503-kinase universe by human UniProt accession or gene symbol). Per condition, every site is robustly
          standardised over the whole measured phosphoproteome: z_i = (x_i - median) / (1.4826 * MAD). Kinase score =
          Stouffer z = sum(z_i) / sqrt(n) over its measured substrates; reported only if n >= 5. Also mean raw log2FC.
Outputs (results/infra_compendium/):
  activity_long.parquet            condition, kinase, method, z, n_substrates, [mean_log2fc, z_core, z_sig, mean, sd]
  activity_corev3_z_wide.parquet   kinase x condition (z), and *_nsub_wide.parquet (n_substrates)
  activity_pspz_z_wide.parquet     kinase x condition (z), and *_nsub_wide.parquet
  psp_ks_edges_human_mapped.tsv    the PSP-only library used for method B (kinase, site)
  activity_coverage.tsv            per kinase: n conditions scored by each method
"""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from common import CACHE, DATA, RES, kinase_universe, read_psp

U = kinase_universe()
by_acc = dict(zip(U["uniprot"], U["kinase"]))
by_name = {k.upper(): k for k in U["kinase"]}

# ---------------------------------------------------------------- PSP-only library
ks = read_psp(DATA / "Kinase_Substrate_Dataset")
orgs = {"human", "mouse", "rat"}
ks = ks[ks["KIN_ORGANISM"].isin(orgs) & ks["SUB_ORGANISM"].isin(orgs)].copy()
k_acc = ks["KIN_ACC_ID"].str.split("-").str[0].map(by_acc).where(ks["KIN_ORGANISM"].eq("human"))
ks["kinase"] = k_acc.fillna(ks["GENE"].str.upper().map(by_name))
ks = ks[ks["kinase"].notna()]
rsd = ks["SUB_MOD_RSD"].str.extract(r"^([STY]\d+)")[0]
human = ks["SUB_ORGANISM"].eq("human") & ~ks["SUB_ACC_ID"].str.contains("-") & rsd.notna()
ks["site"] = (ks["SUB_ACC_ID"] + "_" + rsd).where(human)
ps = read_psp(DATA / "Phosphorylation_site_dataset")
ps = ps[ps["ORGANISM"].eq("human") & ~ps["ACC_ID"].str.contains("-")]
r2 = ps["MOD_RSD"].str.extract(r"^([STY]\d+)")[0]
ps = ps.assign(site=ps["ACC_ID"] + "_" + r2)[r2.notna()]
grp = ps.drop_duplicates("SITE_GRP_ID").set_index("SITE_GRP_ID")["site"]
ks["site"] = ks["site"].fillna(ks["SITE_GRP_ID"].map(grp))
lib = ks.dropna(subset=["site"])[["kinase", "site"]].drop_duplicates()
lib.to_csv(RES / "psp_ks_edges_human_mapped.tsv", sep="\t", index=False)
print(f"PSP-only library: {len(lib)} edges, {lib.kinase.nunique()} kinases, {lib.site.nunique()} sites")

# ---------------------------------------------------------------- method B
sch = pq.read_schema(RES / "site_matrix.parquet").names
conds = sch[1:]
M = pq.read_table(RES / "site_matrix.parquet").to_pandas().set_index("site")
X = M.to_numpy(dtype=np.float32)
med = np.nanmedian(X, axis=0)
mad = np.nanmedian(np.abs(X - med), axis=0) * 1.4826
Z = (X - med) / mad
lib = lib[lib["site"].isin(M.index)]
sites = sorted(lib["site"].unique())
kins = sorted(lib["kinase"].unique())
si = pd.Index(M.index).get_indexer(sites)
Zs = Z[si]                                            # psp sites x conds
Xs = X[si]
fin = np.isfinite(Zs)
Kmat = np.zeros((len(kins), len(sites)), dtype=np.float32)
Kmat[pd.Index(kins).get_indexer(lib["kinase"]), pd.Index(sites).get_indexer(lib["site"])] = 1.0
N = Kmat @ fin.astype(np.float32)
S = Kmat @ np.where(fin, Zs, 0).astype(np.float32)
SX = Kmat @ np.where(fin, Xs, 0).astype(np.float32)
with np.errstate(invalid="ignore", divide="ignore"):
    zB = np.where(N >= 5, S / np.sqrt(N), np.nan)
    mB = np.where(N >= 5, SX / N, np.nan)
B = pd.DataFrame({"condition": np.repeat(np.array(conds, dtype=object)[None, :], len(kins), 0).ravel(),
                  "kinase": np.repeat(np.array(kins, dtype=object), len(conds)),
                  "z": zB.ravel(), "n_substrates": N.ravel().astype(np.int32), "mean_log2fc": mB.ravel()})
B = B[np.isfinite(B["z"])].assign(method="psp_z")
del X, Z, M

# ---------------------------------------------------------------- method A
A = pd.read_parquet(CACHE / "corev3_long.parquet")
assert A["condition"].nunique() == len(conds), A["condition"].nunique()
A = A[np.isfinite(A["z"])].assign(method="corev3")
L = pd.concat([A, B], ignore_index=True)
L["n_substrates"] = L["n_substrates"].astype(np.int32)
for c in ("z", "mean", "sd", "z_core", "z_sig", "mean_log2fc"):
    L[c] = L[c].astype(np.float32)
L = L[["condition", "kinase", "method", "z", "n_substrates", "mean_log2fc", "z_core", "z_sig", "mean", "sd"]]
L.to_parquet(RES / "activity_long.parquet", index=False, compression="zstd")


def wide(df, col):
    w = df.pivot(index="kinase", columns="condition", values=col).reindex(columns=conds)
    return w.astype(np.float32 if col == "z" else np.float32)


for meth, tag in (("corev3", "corev3"), ("psp_z", "pspz")):
    d = L[L.method == meth]
    wide(d, "z").to_parquet(RES / f"activity_{tag}_z_wide.parquet")
    wide(d, "n_substrates").to_parquet(RES / f"activity_{tag}_nsub_wide.parquet")

# ---------------------------------------------------------------- coverage
cv = []
for meth, d in L.groupby("method"):
    g = d.groupby("kinase").agg(n_cond=("z", "size"), med_nsub=("n_substrates", "median"),
                                n_cond_nsub_ge1=("n_substrates", lambda x: int((x >= 1).sum())),
                                n_cond_nsub_ge5=("n_substrates", lambda x: int((x >= 5).sum())))
    g.columns = [f"{meth}_{c}" for c in g.columns]
    cv.append(g)
cov = pd.concat(cv, axis=1).reindex(sorted(set(U["kinase"]) | set(L["kinase"])))
cov.index.name = "kinase"
cov.to_csv(RES / "activity_coverage.tsv", sep="\t")
nC = len(conds)
summ = {}
for meth, d in L.groupby("method"):
    g = d.groupby("kinase")["z"].size()
    summ[meth] = dict(kinases_any=int(len(g)), kinases_ge50pct=int((g >= 0.5 * nC).sum()),
                      kinases_direct_nsub_ge1_any=int(d.loc[d.n_substrates >= 1, "kinase"].nunique()),
                      kinases_nsub_ge5_ge50pct=int((d[d.n_substrates >= 5].groupby("kinase").size() >= 0.5 * nC).sum()),
                      median_kinases_per_condition=float(d.groupby("condition").size().median()),
                      median_kinases_nsub_ge5_per_condition=float(d[d.n_substrates >= 5].groupby("condition").size().median()))
print(pd.DataFrame(summ).to_string())
pd.DataFrame(summ).to_csv(RES / "activity_coverage_summary.tsv", sep="\t")
# method agreement
w1 = wide(L[L.method == "corev3"], "z"); w2 = wide(L[L.method == "psp_z"], "z")
common = w1.index.intersection(w2.index)
a = w1.loc[common].to_numpy().ravel(); b = w2.loc[common].to_numpy().ravel()
ok = np.isfinite(a) & np.isfinite(b)
print(f"core-v3 vs psp_z: {ok.sum()} shared kinase-condition pairs, Pearson r = {np.corrcoef(a[ok], b[ok])[0,1]:.3f}, "
      f"Spearman = {pd.Series(a[ok]).corr(pd.Series(b[ok]), method='spearman'):.3f}")
