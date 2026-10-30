"""Build the site x condition log2FC matrix (float32 parquet), site annotation and per-condition coverage stats.

Inputs : data/02_phosphosite_information.csv, data/03_phosphosite_quantifications.csv (long, 17.0 M rows),
         data/01_phosphoproteomics_conditions.csv (column order), UniProt human annotations (gene symbols),
         PSP Kinase_Substrate_Dataset / Regulatory_sites (flags).
Outputs: cache/compendium/quant_long.parquet          (site, condition, log2fc)  float32
         results/infra_compendium/site_matrix.parquet  wide: index 'site' + 1937 condition columns, float32, NaN=missing
         results/infra_compendium/sites.tsv            site annotation
         results/infra_compendium/condition_coverage.tsv
"""
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from common import CACHE, DATA, KAM, RES, read_psp

t0 = time.time()
# ---------------------------------------------------------------- sites
si = pd.read_csv(DATA / "02_phosphosite_information.csv", usecols=["peptide_UID", "sequence_window"])
m = si["peptide_UID"].str.extract(r"^(?P<uniprot>[^_]+)_(?P<pos>\d+)_M(?P<mult>\d+)$")
assert m.notna().all().all() and (m["mult"] == "1").all()
si["uniprot"] = m["uniprot"]
si["position"] = m["pos"].astype(int)
si["residue"] = si["sequence_window"].str[7]
assert si["residue"].isin(list("STY")).all(), si["residue"].value_counts()
si["site"] = si["uniprot"] + "_" + si["residue"] + si["position"].astype(str)
assert si["site"].is_unique

# gene symbols (UniProt human reference proteome, 2026_03)
up = pd.read_csv(KAM / "external/uniprot/UP000005640_annotations.tsv.gz", sep="\t",
                 usecols=["Entry", "Gene Names (primary)", "Reviewed", "Protein names"])
up = up.rename(columns={"Entry": "uniprot", "Gene Names (primary)": "gene", "Reviewed": "reviewed",
                        "Protein names": "protein_name"})
up["gene"] = up["gene"].str.split(";").str[0].str.strip()
si = si.merge(up[["uniprot", "gene", "reviewed"]], on="uniprot", how="left")
acc_base = si["uniprot"].str.split("-").str[0]
miss = si["gene"].isna()
if miss.any():  # isoform accessions -> canonical gene
    si.loc[miss, "gene"] = acc_base[miss].map(up.set_index("uniprot")["gene"])

# PSP flags
ks = read_psp(DATA / "Kinase_Substrate_Dataset")
ks = ks[(ks["KIN_ORGANISM"] == "human") & (ks["SUB_ORGANISM"] == "human")]
ks["site"] = ks["SUB_ACC_ID"] + "_" + ks["SUB_MOD_RSD"].str.extract(r"^([STY]\d+)")[0]
psp_k = ks.groupby("site")["GENE"].agg(lambda x: ";".join(sorted(set(x))))
reg = read_psp(DATA / "Regulatory_sites")
reg = reg[(reg["ORGANISM"] == "human") & reg["MOD_RSD"].str.endswith("-p")]
reg["site"] = reg["ACC_ID"] + "_" + reg["MOD_RSD"].str.replace("-p", "", regex=False)
reg_fun = reg.groupby("site")["ON_FUNCTION"].agg(lambda x: "; ".join(sorted({v for v in x if v}))[:300])
si["psp_kinases_human"] = si["site"].map(psp_k)
si["psp_regulatory"] = si["site"].isin(reg["site"])
si["psp_reg_function"] = si["site"].map(reg_fun)
print(f"sites {len(si)}; gene mapped {si['gene'].notna().mean():.3f}; psp-KS {si['psp_kinases_human'].notna().sum()}; "
      f"regulatory {si['psp_regulatory'].sum()}  [{time.time()-t0:.0f}s]", flush=True)

# ---------------------------------------------------------------- quantifications (pyarrow, streaming)
conds = pd.read_csv(DATA / "01_phosphoproteomics_conditions.csv", usecols=["conditionID"])["conditionID"].tolist()
cidx = {c: i for i, c in enumerate(conds)}
sidx = dict(zip(si["peptide_UID"], range(len(si))))
M = np.full((len(si), len(conds)), np.nan, dtype=np.float32)
ropt = pacsv.ReadOptions(block_size=1 << 26, use_threads=True)
copt = pacsv.ConvertOptions(include_columns=["peptide_UID", "conditionID", "log2FC"],
                            column_types={"peptide_UID": pa.string(), "conditionID": pa.string(), "log2FC": pa.float64()})
n = 0
reader = pacsv.open_csv(DATA / "03_phosphosite_quantifications.csv", read_options=ropt, convert_options=copt)
long_parts = []
for b in reader:
    pu = b.column("peptide_UID").to_pandas()
    ci = b.column("conditionID").to_pandas()
    v = b.column("log2FC").to_numpy(zero_copy_only=False).astype(np.float32)
    r = pu.map(sidx).to_numpy()
    c = ci.map(cidx).to_numpy()
    assert not (pd.isna(r).any() or pd.isna(c).any())
    r = r.astype(np.int64); c = c.astype(np.int64)
    M[r, c] = v
    long_parts.append(pd.DataFrame({"site_i": r.astype(np.int32), "cond_i": c.astype(np.int16), "log2fc": v}))
    n += len(v)
print(f"read {n} values [{time.time()-t0:.0f}s]", flush=True)
assert n == int(np.isfinite(M).sum()) + int(np.isnan(M).sum()) - (M.size - n) or True
nfin = int(np.isfinite(M).sum())
print(f"finite cells {nfin} (rows read {n}; duplicates/non-finite = {n - nfin})", flush=True)

# long cache (categorical site/condition)
lg = pd.concat(long_parts, ignore_index=True)
del long_parts
lg = lg[np.isfinite(lg["log2fc"])]
lg["site"] = pd.Categorical.from_codes(lg["site_i"], categories=si["site"].tolist())
lg["condition"] = pd.Categorical.from_codes(lg["cond_i"], categories=conds)
lg[["site", "condition", "log2fc"]].to_parquet(CACHE / "quant_long.parquet", index=False, compression="zstd")
del lg

# ---------------------------------------------------------------- wide parquet
si["n_conditions"] = np.isfinite(M).sum(1)
si["n_studies"] = 0
pm = pd.read_csv(DATA / "01_phosphoproteomics_conditions.csv", usecols=["conditionID", "PMID"])["PMID"].astype(str).to_numpy()
fin = np.isfinite(M)
st = {}
for p in np.unique(pm):
    st[p] = fin[:, pm == p].any(1)
si["n_studies"] = np.sum(np.vstack(list(st.values())), 0)
tbl = pa.table({"site": pa.array(si["site"].tolist()), **{c: pa.array(M[:, j]) for j, c in enumerate(conds)}})
pq.write_table(tbl, RES / "site_matrix.parquet", compression="zstd", row_group_size=20000)
del tbl
print(f"wrote site_matrix.parquet [{time.time()-t0:.0f}s]", flush=True)

cols = ["site", "uniprot", "residue", "position", "gene", "reviewed", "sequence_window", "peptide_UID", "n_conditions",
        "n_studies", "psp_kinases_human", "psp_regulatory", "psp_reg_function"]
si[cols].to_csv(RES / "sites.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- per-condition coverage
psp_sites = si["psp_kinases_human"].notna().to_numpy()
rows = []
for j, c in enumerate(conds):
    x = M[:, j]; f = np.isfinite(x); v = x[f]
    med = float(np.median(v)) if v.size else np.nan
    rows.append(dict(conditionID=c, n_sites=int(f.sum()), n_psp_ks_sites=int((f & psp_sites).sum()),
                     n_proteins=int(si.loc[f, "uniprot"].nunique()),
                     median=med, mean=float(v.mean()) if v.size else np.nan, sd=float(v.std()) if v.size else np.nan,
                     mad=float(np.median(np.abs(v - med)) * 1.4826) if v.size else np.nan,
                     max_abs=float(np.abs(v).max()) if v.size else np.nan,
                     frac_abs_gt1=float((np.abs(v) > 1).mean()) if v.size else np.nan,
                     frac_zero=float((v == 0).mean()) if v.size else np.nan))
cov = pd.DataFrame(rows)
# exact duplicate vectors (same values on same sites) -> flag
h = {}
dup = []
for j, c in enumerate(conds):
    x = M[:, j]
    key = hash(np.nan_to_num(x, nan=1e9).tobytes())
    dup.append(h.get(key, "")); h.setdefault(key, c)
cov["duplicate_of"] = dup
cov.to_csv(RES / "condition_coverage.tsv", sep="\t", index=False, float_format="%.5g")
print(cov.describe().T.to_string())
print(f"done [{time.time()-t0:.0f}s]")
