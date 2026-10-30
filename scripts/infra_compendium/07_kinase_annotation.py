"""Kinase annotation table + predicted-kinase columns for sites.

Outputs:
  results/infra_compendium/kinase_annotation.tsv
  results/infra_compendium/sites.tsv (adds pred_top1_kinase, pred_top1_assignment_prob, pred_top3_kinases)
Sources: core-v3 universe (mvp/model/kinases.tsv: Manning/KinHub group+family), kinaseSubMain results/labels/kinase_list.tsv
(subfamily, pseudokinase flag, residue class), IDG understudied kinase list (github druggablegenome/IDGTargets,
IDG_TargetList_CurrentVersion.json, downloaded 2026-10-29 to cache/compendium/external/idg/), PSP (counts),
core-v3 edges, kinaseSubMain release predictions (kinase_top200 / site_top5).
"""
import json

import numpy as np
import pandas as pd

from common import CACHE, KSM, RES, SCRIPTS, kinase_universe

U = kinase_universe()
kl = pd.read_csv(KSM / "results/labels/kinase_list.tsv", sep="\t")
kl = kl[kl["include"]].rename(columns={"gene": "kinase"})
T = U.merge(kl[["kinase", "manning_subfamily", "is_pseudokinase", "residue_class", "protein_name", "hgnc_id"]],
            on="kinase", how="left")
IDG_RENAME = {"ICK": "CILK1", "PRPF4B": "PRP4K", "STK19": "WHR1"}  # current HGNC symbols used by the universe
idg = {IDG_RENAME.get(x["Gene"], x["Gene"]) for x in json.load(open(CACHE / "external/idg/IDG_TargetList_CurrentVersion.json"))
       if x["IDGFamily"] == "Kinase"}
T["idg_dark_kinase"] = T["kinase"].isin(idg)
missing_idg = sorted(idg - set(T["kinase"]))

sites = pd.read_csv(RES / "sites.tsv", sep="\t")
comp = set(sites["site"])
psp = pd.read_csv(RES / "psp_ks_edges_human_mapped.tsv", sep="\t")
T["n_psp_substrate_sites"] = T["kinase"].map(psp.groupby("kinase").size()).fillna(0).astype(int)
T["n_psp_substrate_sites_in_compendium"] = T["kinase"].map(psp[psp.site.isin(comp)].groupby("kinase").size()).fillna(0).astype(int)
e = pd.read_csv(SCRIPTS / "kinact_mvp_pkg/model/edges.tsv.gz", sep="\t")
T["n_corev3_curated_edges"] = T["kinase"].map(e[e.source == "cur"].groupby("kinase").size()).fillna(0).astype(int)
T["n_corev3_curated_edges_in_compendium"] = T["kinase"].map(e[(e.source == "cur") & e.site.isin(comp)].groupby("kinase").size()).fillna(0).astype(int)
T["n_corev3_predicted_edges_in_compendium"] = T["kinase"].map(e[(e.source != "cur") & e.site.isin(comp)].groupby("kinase").size()).fillna(0).astype(int)
top = pd.read_csv(KSM / "results/release/kinase_top200.tsv.gz", sep="\t")
T["n_ksm_top200_sites_in_compendium"] = T["kinase"].map(top[top.site_id.isin(comp)].groupby("gene").size()).fillna(0).astype(int)
cov = pd.read_csv(RES / "activity_coverage.tsv", sep="\t").set_index("kinase")
for c in ["corev3_n_cond", "corev3_n_cond_nsub_ge1", "corev3_n_cond_nsub_ge5", "psp_z_n_cond"]:
    T[c] = T["kinase"].map(cov[c]).fillna(0).astype(int)
T["understudied_by_psp"] = T["n_psp_substrate_sites"] < 5
T = T.rename(columns={"uniprot": "uniprot"})
T.to_csv(RES / "kinase_annotation.tsv", sep="\t", index=False)
print(T["manning_group"].value_counts().to_dict())
print(f"kinases {len(T)}; IDG dark in universe {T.idg_dark_kinase.sum()} (IDG kinases not in universe: {len(missing_idg)}: {missing_idg})")
print(f"PSP substrates >=5: {(T.n_psp_substrate_sites>=5).sum()}, >=5 in compendium: {(T.n_psp_substrate_sites_in_compendium>=5).sum()}")
print(f"dark kinases with psp_z score in >=1 condition: {((T.idg_dark_kinase)&(T.psp_z_n_cond>0)).sum()}; "
      f"with core-v3 n_sub>=1 in >=50% conditions: {((T.idg_dark_kinase)&(T.corev3_n_cond_nsub_ge1>=0.5*1937)).sum()}")

# predicted kinases per compendium site (kinaseSubMain release site_top5)
s5 = pd.read_csv(KSM / "results/release/site_top5.tsv.gz", sep="\t", usecols=["site_id", "rank", "gene", "assignment_prob"])
s5 = s5[s5.site_id.isin(comp) & s5["rank"].le(3)]
t1 = s5[s5["rank"] == 1].set_index("site_id")
t3 = s5.sort_values(["site_id", "rank"]).groupby("site_id")["gene"].agg(";".join)
sites["pred_top1_kinase"] = sites["site"].map(t1["gene"])
sites["pred_top1_assignment_prob"] = sites["site"].map(t1["assignment_prob"]).round(4)
sites["pred_top3_kinases"] = sites["site"].map(t3)
sites.to_csv(RES / "sites.tsv", sep="\t", index=False)
print(f"sites with predictions: {sites.pred_top1_kinase.notna().sum()} / {len(sites)}")
