"""I28 step 1: per-experiment responsiveness of confident kinase-family substrates under verified acute inhibition.
Inputs: I26 LoF conditions + robust z (cache/b_I26), I26 gold sets (PSP in-vivo, LT>=3), I03 kinase conditions.
Experiment ('unit') = study x cell line. Conditions: single-agent inhibitors, time <=120 min or unknown.
Verified = control Stouffer over the family's gold substrates <= -2.5 (>=3 sites); for gold sites the control is
recomputed leaving the tested site out (verification never uses the tested site).
ERK: MEK/RAF inhibitor conditions are added (pathway proxy), verified on the ERK gold set.
Substrates: PSP edges (any member kinase) + kinaseSub release predictions (rank-1 kinase in family, assignment_prob>=0.3,
Cantley percentile >=90 for that kinase). Random controls: 3 per substrate, same residue class (ST/Y), matched on number
of verified units measured, drawn from sites not assigned to the family.
Out: cache/b_I28/L_units.parquet (group, site, unit, study, z, ncond, set), conds_used.tsv
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/b_I26"); sys.path.insert(0, os.path.dirname(__file__))
from lib26 import gold_sets
from lib28 import *
sys.path.insert(0, "/workspace/scripts/infra_compendium"); import load

rng = np.random.default_rng(28)
C = load.conditions().set_index("conditionID")
E = pd.read_parquet(f"{C26}/psp_edges_annot.parquet")
Z = pd.read_parquet(f"{C26}/Z_lof.parquet").set_index("site")
GC = pd.read_csv(f"{C26}/lof_conditions.tsv", sep="\t", dtype={"study": str})
G = pd.read_csv(f"{C26}/groups.tsv", sep="\t"); members = G.groupby("group").kinase.apply(list).to_dict()
GC = GC[GC.kind.eq("inhibitor")]
# ERK pathway proxy
erk_extra = GC[GC.group.isin(["MEK", "RAF"])].assign(group="ERK")
GC = pd.concat([GC.assign(proxy=False), erk_extra.assign(proxy=True)]).drop_duplicates(["group", "conditionID"])
GC["time"] = GC.conditionID.map(C.time_min); GC["cell"] = GC.conditionID.map(C.cell_line_norm).fillna("NA")
GC = GC[GC.time.isna() | (GC.time <= 120)]
GC["unit"] = GC.study + "|" + GC.cell
gold = gold_sets(E)

# ---- substrate sets
Eg = E[E.group.notna() & E.in_compendium]
psp = Eg.groupby(["group", "site"]).agg(psp_LT=("LT_LIT", "max"), psp_vivo=("in_vivo", "max")).reset_index()
mot = pd.read_parquet("/workspace/results/b_I05/cache/motif_pct.parquet")
rel = pd.read_csv(f"{KSM}/results/release/site_top5.tsv.gz", sep="\t", usecols=["site_id", "rank", "gene", "assignment_prob"])
rel = rel[(rel["rank"] == 1) & (rel.assignment_prob >= 0.3)]
k2g = dict(zip(G.kinase, G.group))
rel["group"] = rel.gene.map(k2g); rel = rel[rel.group.notna() & ~rel.group.isin(["MEK", "RAF"]) & rel.site_id.isin(mot.index)]
rel["mpct"] = [mot.at[s, k] if k in mot.columns else np.nan for s, k in zip(rel.site_id, rel.gene)]
pred = rel[rel.mpct >= 90].rename(columns={"site_id": "site"})[["group", "site", "assignment_prob", "mpct"]]
S = psp.merge(pred, on=["group", "site"], how="outer")
S["source"] = np.select([S.psp_LT.notna() & S.assignment_prob.notna(), S.psp_LT.notna()], ["psp+pred", "psp"], "pred")
S["gold"] = [s in gold.get(g, ()) for g, s in zip(S.group, S.site)]
# motif percentile = max over family members
fam_mot = {g: [k for k in ks if k in mot.columns] for g, ks in members.items()}
S["motif_pct"] = [mot.loc[s, fam_mot[g]].max() if (s in mot.index and fam_mot[g]) else np.nan for g, s in zip(S.group, S.site)]
S.to_csv(f"{C28}/substrates.tsv", sep="\t", index=False)
print(S.groupby("group").source.value_counts().unstack(fill_value=0).to_string())

# ---- per group: verified conditions, unit means for all sites (vectorised), LOO for gold
out, used = [], []
Zst = Z.index.str.split("_").str[1].str[0].isin(["S", "T"])
for g, d in GC.groupby("group"):
    cs = [c for c in d.conditionID.unique() if c in Z.columns]
    gs = sorted(gold.get(g, set()) & set(Z.index))
    zg = Z.loc[gs, cs]
    n0 = zg.notna().sum(); s0 = zg.sum()
    stf = s0 / np.sqrt(n0.where(n0 > 0))
    worked = (n0 >= 3) & (stf <= -2.5)
    dd = d.drop_duplicates("conditionID").set_index("conditionID").loc[cs]
    used.append(dd.assign(group=g, ctrl_n=n0, ctrl_stouffer=stf, worked=worked).reset_index())
    wc = [c for c in cs if worked[c]]
    sub = S[S.group == g]
    if not wc and len(sub):
        print(g, "no verified conditions"); continue
    units = dd.loc[wc].unit
    # all-site unit means (non-gold)
    Zw = Z[wc]
    M = pd.DataFrame({u: Zw[list(units.index[units == u])].mean(axis=1) for u in units.unique()})
    Nn = pd.DataFrame({u: Zw[list(units.index[units == u])].notna().sum(axis=1) for u in units.unique()})
    # LOO for gold sites: recompute worked per site
    if gs:
        v = zg.values; meas = ~np.isnan(v)
        n1 = n0.values[None, :] - meas; s1 = s0.values[None, :] - np.where(meas, v, 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            ok = (n1 >= 3) & (s1 / np.sqrt(n1) <= -2.5) & meas & worked.values[None, :]  # LOO-verified AND in verified set
        U = dd.unit.values
        for i, s in enumerate(gs):
            for u in units.unique():
                m = ok[i] & (U == u)
                M.at[s, u] = v[i, m].mean() if m.any() else np.nan
                Nn.at[s, u] = m.sum()
    M = M.where(Nn > 0)
    st_of = {u: u.split("|")[0] for u in M.columns}
    def long(sites, tag):
        sites = [s for s in sites if s in M.index]
        x = M.loc[sites].stack().dropna().rename("z").reset_index(); x.columns = ["site", "unit", "z"]
        nn = Nn.loc[sites].stack(); x["ncond"] = nn.reindex(pd.MultiIndex.from_arrays([x.site, x.unit])).values
        return x.assign(group=g, study=x.unit.map(st_of), set=tag)
    out.append(long(list(sub.site), "substrate"))
    # random matched controls
    nu = M.notna().sum(axis=1); nu = nu[nu > 0]
    pool = nu.drop(index=[s for s in S[S.group == g].site if s in nu.index])
    pst = pool.index.str.split("_").str[1].str[0].isin(["S", "T"])
    rnd = []
    for s in sub.site:
        if s not in nu.index:
            continue
        isst = s.split("_")[1][0] in "ST"; k = min(nu[s], 6)
        cand = pool.index[(pst == isst) & (np.minimum(pool.values, 6) == k)]
        if len(cand):
            rnd += list(rng.choice(cand, min(3, len(cand)), replace=False))
    rnd = sorted(set(rnd))
    out.append(long(rnd, "random"))
    print(g, "verified conds", len(wc), "units", M.shape[1], "studies", len(set(st_of.values())), "substrates measured", nu.reindex(sub.site).notna().sum(), flush=True)

L = pd.concat(out, ignore_index=True)
L.to_parquet(f"{C28}/L_units.parquet")
pd.concat(used).to_csv(f"{C28}/conds_used.tsv", sep="\t", index=False)
print(L.groupby(["group", "set"]).site.nunique().unstack().to_string())
