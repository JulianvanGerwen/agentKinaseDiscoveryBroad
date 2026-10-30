"""Task 1c: is 'AKT/CDK/PKC/mTOR-annotated -> MAPK-axis' beyond chance?
Null A (matched random sites): sites that are assay-'unsupported' for the same group but NOT PSP substrates of it, matched on
  residue, group-motif bin and number of measuring studies; rate of a qualifying MAPK-axis alternative.
Null B (permuted inhibitor labels): each group's verified conditions replaced by an equal number drawn at random from the pool of
  all verified LoF conditions (study labels kept); re-run alternative search for the 270 real unsupported units.
Uses I26 caches read-only (cache/b_I26)."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/b_I26"); from lib26 import classify  # read-only reuse of the classifier rule
OUT = "/workspace/results/audit_F/tables"
rng = np.random.default_rng(7)
C26 = "/workspace/cache/b_I26"
MAPKAX = {"MEK", "RAF", "ERK", "ALK/ROS1", "EGFR/ERBB2", "MET", "ABL/KIT/PDGFR"}
G = pd.read_csv(f"{C26}/groups.tsv", sep="\t"); members = G.groupby("group").kinase.apply(list).to_dict()
DOWN = {"MEK": ["MAPK1", "MAPK3", "RPS6KA1", "RPS6KA3", "RPS6KA5", "MKNK1"], "RAF": ["MAPK1", "MAPK3", "RPS6KA1", "RPS6KA3", "RPS6KA5", "MKNK1"],
        "EGFR/ERBB2": ["MAPK1", "MAPK3", "RPS6KA1", "RPS6KA3", "AKT1", "RPS6KB1"], "MET": ["MAPK1", "MAPK3", "RPS6KA1", "AKT1"],
        "ALK/ROS1": ["MAPK1", "MAPK3", "RPS6KA1", "AKT1"], "ABL/KIT/PDGFR": ["MAPK1", "MAPK3", "AKT1"], "MTOR": ["RPS6KB1", "RPS6KB2", "SGK1"],
        "AKT": ["SGK1"], "ERK": ["RPS6KA1", "RPS6KA3", "RPS6KA5", "MKNK1"]}
mot = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet").astype(np.float32)
gmot = {g: mot[[k for k in set(members[g]) | set(DOWN.get(g, [])) if k in mot.columns]].max(axis=1) for g in members}
GM = pd.DataFrame(gmot)  # site x group: best motif pct among members+downstream

L = pd.read_parquet(f"{C26}/L_main.parquet"); L["group"] = L.group.astype(str); L["study"] = L.study.astype(str)
A = classify(L)
sup = A[A.cls == "supported"]
nsup = sup.groupby("site").size()
U = pd.read_csv(f"{C26}/units_testable.tsv", sep="\t")
E = pd.read_parquet(f"{C26}/psp_edges_annot.parquet")
psp_sites = E[E.group.notna()].groupby("group").site.apply(set).to_dict()


def qualifies(site, g0, Asup, nsup_, motif=True, mapk_only=True):
    a = Asup[(Asup.site == site) & (Asup.group != g0) & (Asup.mean_z <= -1.5)]
    if nsup_.get(site, 0) > 2 or a.empty:
        return False
    if mapk_only:
        a = a[a.group.isin(MAPKAX)]
    if motif:
        a = a[[GM.at[site, g] >= 85 if site in GM.index else False for g in a.group]]
    return len(a) > 0


sup_by_site = {s: d for s, d in sup.groupby("site")}
def q(site, g0, motif=True, mapk_only=True):
    d = sup_by_site.get(site)
    if d is None or nsup.get(site, 0) > 2:
        return False
    d = d[(d.group != g0) & (d.mean_z <= -1.5)]
    if mapk_only: d = d[d.group.isin(MAPKAX)]
    if motif: d = d[[(GM.at[site, g] >= 85) if site in GM.index else False for g in d.group]]
    return len(d) > 0

# ---------- Null A
nst = L.groupby("site").study.nunique()
Zall = pd.read_parquet(f"{C26}/Z_lof.parquet").set_index("site")
vol = ((Zall.abs() > 2).sum(axis=1) / Zall.notna().sum(axis=1)).astype(float)
vol = vol[Zall.notna().sum(axis=1) >= 20]
VEDGES = np.quantile(vol.values, [0.25, 0.5, 0.75, 0.9])
print('vol edges', VEDGES)
res = []
for g0 in ["AKT", "CDK1/2/5/9", "PKC", "MTOR"]:
    obs = U[(U.group == g0) & (U.cls == "unsupported")]
    if obs.empty:
        continue
    for VOL in (False, True):
        pool = A[(A.group == g0) & (A.cls == "unsupported") & ~A.site.isin(psp_sites.get(g0, set()))].copy()
        gm0 = GM[g0] if g0 in GM else None
        def strata(df):
            res_ = df.site.str.split("_").str[1].str[0]
            mb = pd.cut(gm0.reindex(df.site).values, [-1, 50, 85, 95, 101], labels=False)
            nb = pd.cut(nst.reindex(df.site).values, [0, 5, 10, 20, 40, 1000], labels=False)
            k = res_.astype(str).values + "_" + pd.Series(mb).astype(str).values + "_" + pd.Series(nb).astype(str).values
            if VOL:
                k = k + "_" + pd.Series(np.digitize(vol.reindex(df.site).fillna(-1).values, VEDGES)).astype(str).values
            return k
        obs = obs.assign(st=strata(obs)); pool = pool.assign(st=strata(pool))
        o_m = np.mean([q(s, g0) for s in obs.site]); o_any = np.mean([q(s, g0, motif=False, mapk_only=False) for s in obs.site])
        o_mn = np.mean([q(s, g0, motif=False) for s in obs.site])
        pool["q_m"] = [q(s, g0) for s in pool.site]; pool["q_mn"] = [q(s, g0, motif=False) for s in pool.site]
        pool["q_any"] = [q(s, g0, motif=False, mapk_only=False) for s in pool.site]
        ps = pool.groupby("st")
        draws = {"q_m": [], "q_mn": [], "q_any": []}
        cnt = obs.st.value_counts()
        if not any(k in ps.groups for k in cnt.index):
            continue
        for _ in range(500):
            smp = pd.concat([ps.get_group(k).sample(n, replace=True, random_state=int(rng.integers(1e9))) for k, n in cnt.items() if k in ps.groups])
            for c in draws:
                draws[c].append(smp[c].mean())
        for c, o in [("q_m", o_m), ("q_mn", o_mn), ("q_any", o_any)]:
            d = np.array(draws[c])
            res.append(dict(null="A_matched_random_sites" + ("+volatility" if VOL else ""), psp_group=g0, statistic={"q_m": "MAPK-axis alt + motif>=85", "q_mn": "MAPK-axis alt (no motif)",
                            "q_any": "any alt group"}[c], n_obs=len(obs), obs_rate=o, null_mean=d.mean(), null_q05=np.quantile(d, .05),
                            null_q95=np.quantile(d, .95), p_ge=(1 + (d >= o).sum()) / (1 + len(d)), n_pool=len(pool),
                            strata_covered=f"{cnt[cnt.index.isin(ps.groups.keys())].sum()}/{len(obs)}"))
        print(pd.DataFrame(res[-3:]).round(3).to_string())

# ---------- Null B: permuted inhibitor labels
Z = pd.read_parquet(f"{C26}/Z_lof.parquet").set_index("site")
W = pd.read_csv("/workspace/results/b_I26/tables/worked_conditions.tsv", sep="\t", dtype={"study": str})
Wv = W[W.worked]
uns = U[U.cls == "unsupported"][["group", "site"]]
S = sorted(set(uns.site) & set(Z.index))
Zs = Z.loc[S]
cond_study = Wv.drop_duplicates("conditionID").set_index("conditionID").study
pool_c = cond_study.index.values


def run(assign):
    """assign: dict group -> list of conditionIDs. returns per (group,site) classification via classify()."""
    rows = []
    for g, cs in assign.items():
        cs = list(cs)
        st = cond_study.reindex(cs).values
        sub = Zs[cs]
        for s_ in np.unique(st):
            m = sub.loc[:, st == s_].mean(axis=1)
            rows.append(pd.DataFrame({"site": m.index, "z": m.values, "group": g, "study": s_}).dropna())
    LL = pd.concat(rows, ignore_index=True)
    return classify(LL)


def stats(AA):
    sp = AA[AA.cls == "supported"]; ns = sp.groupby("site").size(); sb = {s: d for s, d in sp.groupby("site")}
    n_any = n_mapk = n_mapk_m = 0
    for r in uns.itertuples():
        d = sb.get(r.site)
        if d is None or ns.get(r.site, 0) > 2:
            continue
        d = d[(d.group != r.group) & (d.mean_z <= -1.5)]
        if d.empty: continue
        n_any += 1
        dm = d[d.group.isin(MAPKAX)]
        if len(dm): n_mapk += 1
        dmm = dm[[(GM.at[r.site, g] >= 85) if r.site in GM.index else False for g in dm.group]]
        if len(dmm): n_mapk_m += 1
    return n_any, n_mapk, n_mapk_m

real = {g: d.conditionID.tolist() for g, d in Wv.groupby("group")}
obsB = stats(run(real))
print("observed (no LOO):", obsB)
perm = []
for i in range(300):
    assign = {g: rng.choice(pool_c, size=len(cs), replace=False).tolist() for g, cs in real.items()}
    perm.append(stats(run(assign)))
perm = np.array(perm)
for j, nm in enumerate(["any alternative", "MAPK-axis alternative", "MAPK-axis alternative + motif>=85"]):
    res.append(dict(null="B_permuted_labels", psp_group="all unsupported units (270)", statistic=nm, n_obs=len(uns), obs_rate=obsB[j] / len(uns),
                    null_mean=perm[:, j].mean() / len(uns), null_q05=np.quantile(perm[:, j], .05) / len(uns), null_q95=np.quantile(perm[:, j], .95) / len(uns),
                    p_ge=(1 + (perm[:, j] >= obsB[j]).sum()) / (1 + len(perm))))
R = pd.DataFrame(res)
R.to_csv(f"{OUT}/null_tests.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(R.round(3).to_string())
