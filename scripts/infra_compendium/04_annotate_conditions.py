"""Condition annotation table + drug-target table.

Outputs (results/infra_compendium/):
  conditions_annotated.tsv  one row per conditionID (1,937)
  drug_targets.tsv          one row per agent: curated primary targets + Klaeger 2017 + Hijazi 2020 profiles + confidence
  condition_class_counts.tsv
Reused curated resources (READ-ONLY, kinaseActivityMain/results/ground_truth/, built 2026-10-27 by that project's GT
builder): drug_targets_manual.tsv (agent -> primary targets, notes), drug_aliases.tsv, condition_doses.tsv (doses),
stimulus_kinases.tsv (stimulus -> canonical kinases).  Everything else is rule-based here (see RULES below).
"""
import re

import numpy as np
import pandas as pd

from common import DATA, KAM, RES, kinase_universe

GT = KAM / "results/ground_truth"
cond = pd.read_csv(DATA / "01_phosphoproteomics_conditions.csv", dtype=str).drop(columns=["Unnamed: 0"], errors="ignore")
cond = cond.fillna("")
U = kinase_universe()
UNIV = set(U["kinase"])


def norm(s: str) -> str:
    s = s.replace("Ã¶", "o").replace("ö", "o").replace("&mu;", "u").replace("µ", "u").replace("μ", "u")
    return re.sub(r"\s+", " ", s.replace("-", "")).lower()


# ============================================================== drug / agent dictionary
dtm = pd.read_csv(GT / "drug_targets_manual.tsv", sep="\t", dtype=str).fillna("")
ali = pd.read_csv(GT / "drug_aliases.tsv", sep="\t", dtype=str).fillna("")
# agent -> kind for agents NOT in drug_targets_manual (my annotation; class, subclass, note)
NONKIN = {
    "rituximab": ("other", "antibody_CD20", "anti-CD20 antibody (B-cell receptor-like signalling)"),
    "bortezomib": ("stress", "proteotoxic", "proteasome inhibitor"), "carfilzomib": ("stress", "proteotoxic", "proteasome inhibitor"),
    "metformin": ("other", "metabolic_drug", "biguanide; indirect AMPK activation"),
    "neocarzinostatin": ("stress", "DNA_damage", "radiomimetic DSB inducer"), "hydroxyurea": ("stress", "replication_stress", "RNR inhibitor"),
    "doxorubicin": ("stress", "DNA_damage", "topo II poison"), "4NQO": ("stress", "DNA_damage", "UV-mimetic"),
    "carboplatin": ("stress", "DNA_damage", "platinum"), "cisplatin": ("stress", "DNA_damage", "platinum"),
    "temozolomide": ("stress", "DNA_damage", "alkylating"), "bleomycin": ("stress", "DNA_damage", "DSB inducer"),
    "etoposide": ("stress", "DNA_damage", "topo II poison"), "cytarabine": ("stress", "replication_stress", "antimetabolite"),
    "methotrexate": ("stress", "replication_stress", "antifolate"), "palmitic_acid": ("stress", "lipotoxic", "saturated fatty acid"),
    "coumarin": ("other", "small_molecule_other", "natural compound"), "caffeine": ("other", "small_molecule_other", "natural compound (PDE/ATM-ATR weak inh. at mM)"),
    "quercetin": ("other", "small_molecule_other", "flavonoid; promiscuous weak kinase inhibitor"),
    "curcumin": ("other", "small_molecule_other", "natural compound"), "genistein": ("other", "small_molecule_other", "isoflavone; weak broad TK inhibitor"),
    "ginsenoside_rg3": ("other", "small_molecule_other", "natural compound"), "24D": ("other", "small_molecule_other", "2,4-D herbicide"),
    "nicotine": ("other", "receptor_agonist", "nAChR agonist"), "a_bungarotoxin": ("other", "receptor_antagonist", "nAChR antagonist"),
    "romidepsin": ("other", "epigenetic_drug", "HDAC inhibitor"), "vorinostat": ("other", "epigenetic_drug", "HDAC inhibitor"),
    "A485": ("other", "epigenetic_drug", "p300/CBP HAT inhibitor"), "A486": ("other", "epigenetic_drug", "inactive A-485 analogue"),
    "danuglipron": ("other", "receptor_agonist", "GLP-1R agonist (small molecule)"), "semaglutide": ("other", "receptor_agonist", "GLP-1R agonist (peptide)"),
    "iloprost": ("other", "receptor_agonist", "prostacyclin IP-receptor agonist (cAMP/PKA)"),
    "PRL3i_unspecified": ("other", "phosphatase_inhibitor", "PRL-3 phosphatase inhibitor"),
    "pervanadate": ("other", "phosphatase_inhibitor", "pan-PTP inhibitor (tyrosine phosphorylation up)"),
    "okadaic_acid": ("other", "phosphatase_inhibitor", "PP2A/PP1 inhibitor"), "cantharidic_acid": ("other", "phosphatase_inhibitor", "PP2A inhibitor"),
    "SHP099": ("other", "phosphatase_inhibitor", "SHP2 allosteric inhibitor (RAS-ERK down indirectly)"),
    "dyngo": ("other", "small_molecule_other", "dynamin inhibitor (Dyngo-4a)"), "doxycycline": ("genetic", "inducible_expression", "dox-inducible construct"),
    "AZ123": ("other", "small_molecule_other", "unidentified compound (AML study 23532336)"),
    "TRC794": ("other", "small_molecule_other", "unidentified compound TRC-794"),
    "RG7356": ("other", "antibody_other", "anti-CD44 antibody"), "maraviroc": ("other", "receptor_antagonist", "CCR5 antagonist"),
    "noscapine": ("other", "cell_cycle", "microtubule-binding"), "nocodazole": ("other", "cell_cycle", "microtubule depolymeriser (mitotic arrest)"),
    "paclitaxel": ("other", "cell_cycle", "microtubule stabiliser"), "thymidine": ("other", "cell_cycle", "G1/S arrest"),
    "dexmedetomidine": ("other", "receptor_agonist", "alpha2-adrenergic agonist"), "propofol": ("other", "small_molecule_other", "anaesthetic"),
    "tamoxifen": ("other", "small_molecule_other", "SERM"), "olaparib": ("other", "small_molecule_other", "PARP inhibitor"),
    "CB5083": ("other", "small_molecule_other", "p97 inhibitor"), "NMS873": ("other", "small_molecule_other", "p97 inhibitor"),
    "edelfosine": ("other", "small_molecule_other", "alkyl-lysophospholipid (PLC/PKC axis)"), "U73122": ("other", "small_molecule_other", "PLC inhibitor"),
    "cortistatinA": ("kinase_inhibitor", "", "CDK8/19 inhibitor"),
    "MAPKi_unspecified": ("kinase_inhibitor", "", "unspecified 'MAPKi'"),
}
EXTRA_ALIASES = {  # alias -> agent (in addition to drug_aliases.tsv); spellings seen in 01_ metadata
    "Rituximab": "rituximab", "Bortezomib": "bortezomib", "Velcade": "bortezomib", "Carfilzomib": "carfilzomib",
    "metformin": "metformin", "NCS": "neocarzinostatin", "HU": "hydroxyurea", "doxorubicin": "doxorubicin",
    "4-nitroquinoline": "4NQO", "carbo": "carboplatin", "CP24h": "carboplatin", "TMZ": "temozolomide", "Bleomycin": "bleomycin",
    "Bleomcyin": "bleomycin", "etoposide": "etoposide", "Cisplatin": "cisplatin", "Cytarabine": "cytarabine",
    "Methotrexat": "methotrexate", "Palmitic Acid": "palmitic_acid", "coumarin": "coumarin", "caffeine": "caffeine",
    "quercetin": "quercetin", "Curcumin": "curcumin", "Genistein": "genistein", "Ginsenoside": "ginsenoside_rg3", "24D": "24D",
    "nicotine": "nicotine", "A-bungarotoxin": "a_bungarotoxin", "Romidepsin": "romidepsin", "Vorinostat": "vorinostat",
    "A485": "A485", "A486": "A486", "Danuglipron": "danuglipron", "Semaglutide": "semaglutide", "Iloprost": "iloprost",
    "PRL3i": "PRL3i_unspecified", "Pervanadate": "pervanadate", "Dyngo": "dyngo", "Doxycycline": "doxycycline",
    "AZ123": "AZ123", "TRC-794": "TRC794", "RG7356": "RG7356", "maraviroc": "maraviroc", "Noscapine": "noscapine",
    "Nocodazole": "nocodazole", "Taxol": "paclitaxel", "Paclitaxel": "paclitaxel", "Thymidine": "thymidine",
    "Dexmedetomidine": "dexmedetomidine", "Propofol": "propofol", "Tamoxifen": "tamoxifen", "Olaparib": "olaparib",
    "CB-5083": "CB5083", "NMS-873": "NMS873", "Edelfosine": "edelfosine", "U73122": "U73122", "Cortistatin A": "cortistatinA",
    "MAPKi": "MAPKi_unspecified", "SB202190": "SB202190", "Ku-794": "KU0063794", "PI-103": "PI103", "1-NM-PP1": "NM_PP1_1",
    "Plk1(as) 3-MB-PP1": "PLK1as_3MBPP1", "BI 4834": "BI4834", "BI2536": "BI2536", "BI-D1870": "BID1870", "SC-514": "SC514",
    "SM": "SM_unspecified", "LRRK2-IN-1": "LRRK2IN1", "VE-821": "VE821", "KU60019": "KU60019", "ATMi": "ATMi", "ATRi": "ATRi",
    "CHK1i": "CHK1i", "DNA-PKi": "DNAPKi", "DNA-Pki": "DNAPKi", "p38i": "p38i", "BRAFi": "BRAFi", "MEKi": "MEKi",
    "BisindolylmaleimideI": "BIM1", "Okadaic": "okadaic_acid", "Cantharidic": "cantharidic_acid", "SYK inhibition": "SYKi",
    "Akti": None,  # prefix in 32034161 ("Akti_MK-2206"); resolved by compound name
}
NONKIN["SM_unspecified"] = ("other", "small_molecule_other", "'SM' (Smac mimetic?) in 28539327 - unresolved")
GENERIC = set(ali.loc[ali["is_class"].eq("1"), "agent"]) | {"AKTi_unspecified"}
alias = {}
for a, g in zip(ali["alias"], ali["agent"]):
    alias[norm(a)] = g
for a in dtm["agent"]:
    alias[norm(a)] = a
for a, g in EXTRA_ALIASES.items():
    if g is None:
        alias.pop(norm(a), None)
    else:
        alias[norm(a)] = g
for bad in ["bio", "sm", "hu", "ncs"]:   # short aliases: matched only as whole tokens (see find_agents)
    pass
ALIAS_KEYS = sorted(alias, key=len, reverse=True)
SHORT = {k for k in ALIAS_KEYS if len(k) <= 4}


def find_agents(text: str) -> list:
    t = norm(text)
    found = []
    for k in ALIAS_KEYS:
        if k in SHORT:
            pat = r"(?<![a-z0-9])" + re.escape(k) + r"(?![a-z])"
        else:
            pat = r"(?<![a-z])" + re.escape(k) + r"(?![a-z])"
        m = re.search(pat, t)
        if m:
            found.append(alias[k])
            t = t[:m.start()] + "#" * (m.end() - m.start()) + t[m.end():]
    return sorted(set(found))


ZERO_DOSE = re.compile(r"(^|[_\s])0(\.0+)?\s?(nM|uM|μM|ng/mL|mM)\b|^[A-Za-z]+0 μM|_0uM", re.I)
ZERO_TIME = re.compile(r"(^|[_\s])0\s?(h|min|m)(?![A-Za-z])")
CTRL_LIKE = re.compile(r"^(Ctr|control|DMSO|Vehicle|untreated|Mock|WT|Parental)$", re.I)
# PMID-level overrides: pmid -> (class, subclass, stimulus, stim_class)
PMID_OVERRIDE = {"34974192": ("stress", "mito_depolarisation", "mitochondrial depolarisation (parkin-expressing HeLa)", "MITO_DEPOL")}

# ============================================================== stimuli / other perturbation rules (regex on text)
# (regex, class, subclass, stimulus_name, stim_class for stimulus_kinases.tsv)
STIM = [
    (r"\bHB-EGF\b", "growth_factor", "RTK_ligand", "HB-EGF", "EGF"),
    (r"(?<![A-Za-z-])EGF(?![A-Za-z])", "growth_factor", "RTK_ligand", "EGF", "EGF"),
    (r"FGF(7|10|2)", "growth_factor", "RTK_ligand", "FGF", "FGF"),
    (r"\bNGF\b", "growth_factor", "RTK_ligand", "NGF", "NGF"),
    (r"[Ii]nsulin", "growth_factor", "RTK_ligand", "insulin", "INSULIN"),
    (r"[Ss]erum", "growth_factor", "serum", "serum", ""),
    (r"TNF", "cytokine", "TNF", "TNF", "TNF"),
    (r"IFNg/LPS", "cytokine", "IFNG_LPS", "IFNg+LPS", "IFNG_LPS"),
    (r"(?<![A-Za-z])IFN", "cytokine", "interferon", "IFN", "IFN"),
    (r"IL4", "cytokine", "interleukin", "IL-4", "IL4"),
    (r"IL-23", "cytokine", "interleukin", "IL-23", "IL23"),
    (r"TGFB", "cytokine", "TGFb", "TGF-beta", "TGFB"),
    (r"CXCL12|SDF-1", "cytokine", "chemokine", "CXCL12", ""),
    (r"Anti-CD3|CD3\+CD4", "other", "immune_receptor", "TCR stimulation", "TCR"),
    (r"BCR stimulation", "other", "immune_receptor", "BCR stimulation", "BCR"),
    (r"\bLPA\b", "other", "GPCR_ligand", "LPA", ""),
    (r"\bAng\b", "other", "GPCR_ligand", "angiotensin II", ""),
    (r"SII Ang", "other", "GPCR_ligand", "SII-angiotensin (biased)", ""),
    (r"Thrombin", "other", "GPCR_ligand", "thrombin", ""),
    (r"PGE2", "other", "GPCR_ligand", "PGE2", "CAMP_PKA"),
    (r"EDN", "other", "GPCR_ligand", "endothelin", ""),
    (r"\bADP\b", "other", "GPCR_ligand", "ADP", ""),
    (r"\bPMA\b", "other", "PKC_activator", "PMA", "PMA"),
    (r"\bNCM\b", "other", "other_stimulus", "NCM (conditioned medium)", ""),
    (r"\bPOS\b", "other", "other_stimulus", "photoreceptor outer segments", ""),
    (r"xray|proton|carbon|\bIR\b|Irradiated", "stress", "DNA_damage", "ionising radiation", "DSB"),
    (r"UVB?(?![A-Za-z])", "stress", "UV", "UV", "UV"),
    (r"H2O2", "stress", "oxidative", "H2O2", ""),
    (r"\bcold\b", "stress", "cold_ischemia", "cold ischaemia", ""),
    (r"Starv|SFM|NoSerum", "stress", "starvation", "starvation/serum-free", ""),
    (r"\bHG\b", "stress", "metabolic", "high glucose", ""),
    (r"\bAO\b", "stress", "mito_depolarisation", "antimycin/oligomycin", "MITO_DEPOL"),
    (r"[Pp]arkin|PARKIN", "genetic", "overexpression", "parkin expression (+ mito depolarisation)", "MITO_DEPOL"),
    (r"SARS|HIV|HSV|flexneri|ZIKV|IAV|trachomatis|strain infected|aureus|cocult", "other", "infection_or_coculture", "", ""),
    (r"TcdB|tBID", "other", "toxin", "", ""),
    (r"phase|^G1|^G2$|^M$|^S$|Early S|Late S|[Mm]itosis|Mitotic|Asynchronous|G1_S|Cyclin A", "other", "cell_cycle", "cell-cycle state", ""),
    (r"Exercise", "other", "physiology", "exercise", ""),
    (r"3D|2D", "other", "culture_format", "3D culture", ""),
    (r"OIS", "other", "senescence", "oncogene-induced senescence", ""),
]
# genetic perturbations: regex -> (subclass, kinase(s) perturbed in the 503 universe or '', direction of kinase change)
GEN = [
    (r"MARK3_KO", "knockout", "MARK3", -1), (r"EDNRB_KO", "knockout", "", 0), (r"PTPN11_KO", "knockout", "", 0),
    (r"sh-TBK1", "knockdown", "TBK1", -1), (r"PINK1 KD", "knockdown", "PINK1", -1), (r"PINK1 WT", "rescue_vs_KD", "PINK1", 1),
    (r"PKD1kd", "knockdown", "PRKD1", -1), (r"RSK1/2 RNAi", "knockdown", "RPS6KA1;RPS6KA3", -1),
    (r"Fam20C KO", "knockout", "FAM20C", -1), (r"AMPK-DKO", "knockout", "PRKAA1;PRKAA2", -1),
    (r"PKD2\+PKD3 KD", "knockdown", "PRKD2;PRKD3", -1), (r"PKD2 KD", "knockdown", "PRKD2", -1),
    (r"CDKL5 expressing", "rescue_vs_KO", "CDKL5", 1), (r"MASTL overexpression", "overexpression", "MASTL", 1),
    (r"IKKB\(K44M\)", "dominant_negative", "IKBKB", -1), (r"IKKB\(WT\)", "overexpression", "IKBKB", 1),
    (r"BGLF4|Vector", "viral_kinase_vs_vector", "", 0), (r"Raptor KD", "knockdown", "", 0),
    (r"si(CT45|PP4C)", "knockdown", "", 0), (r"CD79a KD", "knockdown", "", 0), (r"CHC siRNA", "knockdown", "", 0),
    (r"K8-S", "phosphosite_mutant", "", 0), (r"Ex9-KI|Ex20-KI", "knock_in", "", 0),
    (r"IDH1mut|RAS-NHA", "oncogene_expression", "", 0), (r"Cas9", "crispr", "", 0), (r"TOPBP1-AID|ETAA1-AAD", "degron_or_mutant", "", 0),
    (r"CD38-trans", "overexpression", "", 0), (r"C431A-parkin|WT-parkin|ParkinS65A", "overexpression", "", 0),
    (r"HKD|EKD", "knockdown", "", 0), (r"Vpr-deleted", "viral_mutant", "", 0),
    (r"^A4$|^H6$", "clone_comparison", "", 0), (r"\bABL\b", "fusion_vs_control", "ABL1", 1), (r"^GLP_", "overexpression", "", 0),
]
DISEASE = (r"HCC|[Cc]ancer|Alzheimer|carcinoma|Metastatic|SW620|JHU|FADU|CAL27|ARMS|ERMS|^G4$|^WNT$|resistant|"
           r"sensitive|DasR|PazR|HER2\+|^AML$|P31/Fuj|^CAF$|Epithelial|VMCUB|SW780|RT112|UMUC3|T24|J82|SCC-9|"
           r"C4-2|DU145|EF-1|LNCaP|NCI-H660|PARCB|PC3|VCaP|22Rv1")

NO_NET = {"infection_or_coculture", "toxin", "cell_cycle", "culture_format", "senescence", "physiology"}
TIME_RE = re.compile(r"(?<![A-Za-z0-9.])(\d+(?:\.\d+)?)\s*(sec|s|min|m|hrs|hr|h|days|day|d)(?![A-Za-z])")
UNIT = {"sec": 1 / 60, "s": 1 / 60, "min": 1, "m": 1, "h": 60, "hr": 60, "hrs": 60, "d": 1440, "day": 1440, "days": 1440}


def parse_time(s: str):
    s = s.strip()
    if re.fullmatch(r"\d+(\.\d+)?", s):          # Ochoa2016 rows: bare numbers are minutes
        return float(s)
    m = TIME_RE.fullmatch(s)
    return float(m.group(1)) * UNIT[m.group(2)] if m else np.nan


def time_tokens(s: str):
    return [float(a) * UNIT[u] for a, u in TIME_RE.findall(s)]


def strip_time(s: str) -> str:
    s = TIME_RE.sub("", s)
    s = re.sub(r"(?<=\d)h(?=ng/mL)", "", s)       # decryptM 'Rituximab_1hng/mL'
    return re.sub(r"[_\s]+", " ", s).strip(" _,")


def cell_norm(s: str) -> str:
    k = re.sub(r"[-\s.]", "", s).upper()
    k = {"MCF7": "MCF7", "HELA013": "HELA", "HELAS3": "HELA", "HELAKYOTO": "HELA", "HELAFUCCI": "HELA",
         "U2OS": "U2OS", "HEPG2": "HEPG2", "K562": "K562", "HEK293": "HEK293", "HEK293T": "HEK293T"}.get(k, k)
    return k


TISSUE = re.compile(r"liver|brain|muscle|breast cancer|pancreas|blood|plasma|medulloblastoma|HNSCC|colorectal cancer|"
                    r"Rhabdomyosarcoma|islet|Bone Marrow|carcinoma|organoid|CD4\+|BOEC|Platelet|cancer stem|T cell", re.I)

# ============================================================== annotate
stim_k = pd.read_csv(GT / "stimulus_kinases.tsv", sep="\t", comment="#")
stim_direct = stim_k[stim_k["cls"] == "direct"].groupby("stim_class")["kinase"].agg(lambda x: ";".join(x))
dose = pd.read_csv(GT / "condition_doses.tsv", sep="\t")
dose_c = dose[dose.side == "condition"].groupby("conditionID")["dose_nM"].max()
dt = dtm.set_index("agent")

rows = []
for r in cond.itertuples(index=False):
    c, ctl, tp, smp = r.Condition, r.Control, r.Timepoint, r.Sample
    # --- time
    tmin = parse_time(tp)
    if np.isnan(tmin) and not np.isnan(parse_time(smp)) and tp:      # swapped Sample/Timepoint (37783352)
        tmin, smp = parse_time(smp), tp
    toks = time_tokens(c)
    tsrc = "Timepoint" if not np.isnan(tmin) else ""
    if len(set(toks)) == 1 and (np.isnan(tmin) or abs(toks[0] - tmin) > 1e-6):
        tmin, tsrc = toks[0], "condition_text"
    tnote = "long-term/steady-state" if tp.strip().lower() == "long" else ""
    # --- agents (net of control), with zero-dose / time-zero controls treated as untreated baseline
    ctl_eff = "" if (ZERO_DOSE.search(ctl) or ZERO_TIME.search(ctl)) else ctl
    ca = find_agents(c)
    ka = find_agents(ctl_eff)
    reversed_ = False
    is_tc = (strip_time(c).lower() == strip_time(ctl).lower() and bool(time_tokens(ctl)) and time_tokens(ctl)[0] > 0) or \
        re.fullmatch(r"\s*\d+\s*(min|h)\s*", ctl or "") is not None or ctl.strip() == "Rituximab_1min"
    net = ca if is_tc else [a for a in ca if a not in ka]
    background = [a for a in ca if a in ka] if not is_tc else []
    extra_ctl = [a for a in ka if a not in ca]
    if not net and extra_ctl and not is_tc:          # e.g. 'DMSO vs Ceritinib', 'Bleomycin DMSO vs Bleomycin DNA-PKi'
        net, reversed_ = extra_ctl, True
    kin_inh = [a for a in net if (a in dt.index and dt.loc[a, "mechanism"] not in ("non_kinase",)) or
               (a in NONKIN and NONKIN[a][0] == "kinase_inhibitor")]
    if any(a not in GENERIC for a in net):           # specific compound named -> drop class-level alias hits
        net = [a for a in net if a not in GENERIC]
        kin_inh = [a for a in kin_inh if a not in GENERIC]
    other_drugs = [a for a in net if a not in kin_inh]
    # --- stimuli (net of control text, except context-type rules)
    stims = []
    for rx, cl, sub, name, sc in STIM:
        if sub in NO_NET:
            hit = re.search(rx, c) or re.search(rx, ctl)
        else:
            hit = re.search(rx, c) and not re.search(rx, ctl_eff)
        if hit:
            stims.append((cl, sub, name or hit.group(0), sc))
    if not stims and CTRL_LIKE.match(c.strip()):           # 'Ctr vs Serum' -> reversed stimulus
        for rx, cl, sub, name, sc in STIM:
            if re.search(rx, ctl_eff):
                stims.append((cl, sub, name or re.search(rx, ctl_eff).group(0), sc)); reversed_ = True
    gens = [(sub, k, d) for rx, sub, k, d in GEN if re.search(rx, c) and not re.search(rx, ctl)]
    gens_rev = [] if gens else [(sub, k, -d) for rx, sub, k, d in GEN if re.search(rx, ctl) and not re.search(rx, c)]
    if not gens and not gens_rev:
        gens = [(sub, "", 0) for rx, sub, k, d in GEN if re.search(rx, c)][:1]
    disease = bool(re.search(DISEASE, c)) or bool(re.search(r"^Normal$|^control$|OKF6|BTERT|GMPB|HaCaT|NPF|Mesenchymal|HER2-|Non-metastatic|insensitive|^Resistant$|subgroup|Kasumi", ctl)
                                                     and re.search(DISEASE, c + " " + smp))
    # --- class (priority: kinase inhibitor > genetic > growth factor/cytokine > stress > disease > other)
    sub, stim_name, stim_cls = "", "", ""
    if kin_inh:
        pclass = "kinase_inhibitor"
        sub = "inhibitor_reversed_contrast" if reversed_ else ("inhibitor+stimulus" if stims else "inhibitor")
        ctl_kin = [a for a in ka if a not in ca and ((a in dt.index and dt.loc[a, "mechanism"] != "non_kinase")
                                                     or NONKIN.get(a, ("",))[0] == "kinase_inhibitor")]
        if ctl_kin and not reversed_:
            sub = "inhibitor_vs_other_inhibitor"
    elif gens or gens_rev:
        pclass = "genetic"; g = (gens or gens_rev)[0]; sub = g[0]
    elif any(s[0] in ("growth_factor", "cytokine") for s in stims):
        s = [s for s in stims if s[0] in ("growth_factor", "cytokine")][0]
        pclass, sub, stim_name, stim_cls = s
    elif other_drugs and NONKIN.get(other_drugs[0], ("other",))[0] == "stress":
        pclass, sub = "stress", NONKIN[other_drugs[0]][1]
    elif stims and stims[0][0] == "stress":
        pclass, sub, stim_name, stim_cls = stims[0]
    elif other_drugs:
        k0 = NONKIN.get(other_drugs[0], ("other", "small_molecule_other", ""))
        pclass, sub = k0[0], k0[1]
    elif stims:
        pclass, sub, stim_name, stim_cls = stims[0]
    elif disease:
        pclass, sub = "disease_vs_control", "tumour_or_cell_line_comparison"
    else:
        pclass, sub = "other", "unclassified"
    if r.PMID in PMID_OVERRIDE:
        pclass, sub, stim_name, stim_cls = PMID_OVERRIDE[r.PMID]
    # same treatment, later vs earlier (non-zero) time -> time_course contrast
    if is_tc:
        pclass, sub = "time_course", f"vs_earlier_time ({sub})"
    # --- targets
    tgt, tconf, tsrc_t = [], "", ""
    if kin_inh:
        for a in kin_inh:
            if a in dt.index:
                tgt += [k for k in dt.loc[a, "direct_targets"].split(";") if k]
    gkin = []
    for g in (gens or gens_rev):
        if g[1]:
            gkin += [f"{k}({'+' if g[2] > 0 else '-'})" for k in g[1].split(";")]
    if stims and not stim_cls:
        stim_cls = next((s[3] for s in stims if s[3]), "")
    if not stim_name and stims:
        stim_name = "+".join(dict.fromkeys(s[2] for s in stims))
    rows.append(dict(
        conditionID=r.conditionID, study_pmid=r.PMID, source=r.source, condition_raw=c, control_raw=ctl,
        perturbation_class=pclass, perturbation_subclass=sub,
        agents=";".join(net), kinase_inhibitors=";".join(kin_inh), other_drugs=";".join(other_drugs),
        control_agents=";".join(ka), background_agents=";".join(background), contrast_reversed=reversed_,
        stimulus=stim_name, stimulus_class=stim_cls,
        stimulus_expected_up_kinases=stim_direct.get(stim_cls, "") if stim_cls else "",
        inhibitor_primary_targets=";".join(dict.fromkeys(tgt)),
        genetic_target_kinases=";".join(gkin),
        dose_nM=dose_c.get(r.conditionID, np.nan) if kin_inh or other_drugs else np.nan,
        cell_type=smp, cell_line_norm=cell_norm(smp),
        sample_type="tissue/primary/organoid" if TISSUE.search(smp) else "cell_line",
        time_raw=tp, time_min=tmin, time_source=tsrc, time_note=tnote,
        labeling=r[cond.columns.get_loc("Labeling.method")], spectrometer=r.Spectrometer))
A = pd.DataFrame(rows)

# ============================================================== drug-target table (Klaeger + Hijazi)
# Klaeger 2017 Table S2 'Kinobeads' (4-cell-line lysate): high-confidence targets with apparent Kd
kl = pd.read_excel(DATA / "Klaeger/aan4368_Table_S2.xlsx", sheet_name="Kinobeads")
kl = kl[kl["Target Classification"].eq("High confidence")]
kl = kl.assign(gene=kl["Gene Name"].astype(str).str.split(";")).explode("gene")
kl_kd = kl.groupby(["Drug", "gene"])["Apparent Kd"].min()
kl_drugs = sorted(kl["Drug"].unique())
kl_des = pd.read_excel(DATA / "Klaeger/aan4368_Table_S1.xlsx", sheet_name="Inhibitor annotation").set_index("Drug")
kl_norm = {norm(d): d for d in kl_drugs}
hz = pd.read_csv(DATA / "kinhib_inhib_comb__1.csv", index_col=0)
hz_norm = {}
for col in hz.columns:
    a = alias.get(norm(col))
    if a:
        hz_norm[a] = col
agents_used = sorted({a for s in A["agents"] for a in s.split(";") if a} | {a for s in A["control_agents"] for a in s.split(";") if a})
drows = []
for a in agents_used:
    rec = dict(agent=a)
    if a in dt.index:
        d = dt.loc[a]
        rec.update(canonical_name=d["canonical_name"], kind="kinase_inhibitor" if d["mechanism"] != "non_kinase" else "other",
                   mechanism=d["mechanism"], primary_targets=d["direct_targets"], curated_note=d["note"])
        kn = d["klaeger_name"] or ""
    else:
        k0 = NONKIN.get(a, ("other", "", ""))
        rec.update(canonical_name=a, kind=k0[0], mechanism=k0[1], primary_targets="", curated_note=k0[2])
        kn = ""
    kcol = kl_norm.get(norm(kn)) or kl_norm.get(norm(a)) or kl_norm.get(norm(rec["canonical_name"]))
    if kcol:
        v = kl_kd.loc[kcol].sort_values() if kcol in kl_kd.index.get_level_values(0) else pd.Series(dtype=float)
        rec["klaeger_drug"] = kcol
        rec["klaeger_designated"] = str(kl_des.loc[kcol, "Designated targets"]) if kcol in kl_des.index else ""
        rec["klaeger_kd_le100nM"] = ";".join(v[v <= 100].index)
        rec["klaeger_kd_le1000nM"] = ";".join(v[v <= 1000].index)
    hcol = hz_norm.get(a)
    if hcol:
        rec["hijazi_compound"] = hcol
        rec["hijazi_targets"] = ";".join(hz.index[hz[hcol] == 1])
    prim = set(filter(None, rec["primary_targets"].split(";")))
    prof = set(filter(None, (rec.get("klaeger_kd_le1000nM", "") + ";" + rec.get("hijazi_targets", "")).split(";")))
    note = (rec["curated_note"] or "").lower()
    if rec["kind"] != "kinase_inhibitor" or not prim:
        conf = "none" if rec["kind"] != "kinase_inhibitor" else "low"
    elif "inferred" in note or "not named" in note or "assumed" in note or rec["mechanism"] in ("indirect", "antibody", "chemical_genetic", "pan"):
        conf = "low"
    elif prim & prof:
        conf = "high"
    elif rec["mechanism"] == "lipid_only":
        conf = "medium"
    else:
        conf = "medium"
    rec["target_confidence"] = conf
    rec["primary_targets_in_universe"] = ";".join(k for k in rec["primary_targets"].split(";") if k in UNIV)
    rec["profiled_targets_in_universe"] = ";".join(sorted(k for k in prof if k in UNIV))
    drows.append(rec)
D = pd.DataFrame(drows)
D.to_csv(RES / "drug_targets.tsv", sep="\t", index=False)
dconf = D.set_index("agent")["target_confidence"]
dprof = D.set_index("agent")["profiled_targets_in_universe"].fillna("")
CONF_RANK = {"high": 3, "medium": 2, "low": 1, "none": 0}


def cond_conf(s):
    ags = [a for a in s.split(";") if a]
    if not ags:
        return ""
    return min((dconf.get(a, "low") for a in ags), key=lambda x: CONF_RANK[x])


A["target_confidence"] = A["kinase_inhibitors"].map(cond_conf)
A.loc[A["perturbation_subclass"].eq("inhibitor_reversed_contrast") & A["target_confidence"].eq("high"), "target_confidence"] = "medium"
A.loc[A["perturbation_subclass"].eq("inhibitor_vs_other_inhibitor"), "target_confidence"] = "low"
A["inhibitor_profiled_targets"] = A["kinase_inhibitors"].map(
    lambda s: ";".join(sorted({k for a in s.split(";") if a for k in dprof.get(a, "").split(";") if k})))
A["inhibitor_primary_targets_in_universe"] = A["inhibitor_primary_targets"].map(
    lambda s: ";".join(k for k in s.split(";") if k in UNIV))
A["is_combination"] = A["agents"].str.count(";").ge(1)

# ============================================================== QC from coverage
cov = pd.read_csv(RES / "condition_coverage.tsv", sep="\t")
A = A.merge(cov[["conditionID", "n_sites", "n_psp_ks_sites", "median", "mad", "duplicate_of"]], on="conditionID", how="left")
A["qc_flag"] = np.select([A["duplicate_of"].fillna("").ne(""), A["mad"].gt(3) | A["median"].abs().gt(1), A["n_sites"].lt(3000)],
                         ["duplicate_vector", "scale_suspect", "low_coverage"], "ok")
A = A.drop(columns=["median", "mad"])

# ============================================================== time-series & dose-series groups
A["treatment_key"] = A["condition_raw"].map(strip_time).str.lower() + " | " + A["control_raw"].map(strip_time).str.lower()
grp = A.groupby(["study_pmid", "cell_line_norm", "treatment_key"])
A["_nt"] = grp["time_min"].transform(lambda x: x.dropna().nunique())
ts = A[A["_nt"] >= 2].groupby(["study_pmid", "cell_line_norm", "treatment_key"]).ngroup()
A["time_series_id"] = ""
A.loc[ts.index, "time_series_id"] = ["TS%04d" % i for i in ts]
A["time_series_n_timepoints"] = np.where(A["time_series_id"].ne(""), A["_nt"], 0)
A["time_rank"] = np.where(A["time_series_id"].ne(""),
                          A.groupby("time_series_id")["time_min"].rank(method="dense"), np.nan)
dkey = A["agents"].where(A["agents"].ne(""))
A["_dk"] = A["study_pmid"] + "|" + A["cell_line_norm"] + "|" + dkey.fillna("") + "|" + A["time_min"].astype(str)
nd = A[dkey.notna()].groupby("_dk")["dose_nM"].transform("nunique")
A["dose_series_id"] = ""
sel = nd[nd >= 3].index
A.loc[sel, "dose_series_id"] = "DS" + A.loc[sel].groupby("_dk").ngroup().map("{:04d}".format)
A = A.drop(columns=["_nt", "_dk"])

order = ["conditionID", "study_pmid", "source", "perturbation_class", "perturbation_subclass", "condition_raw", "control_raw",
         "agents", "kinase_inhibitors", "other_drugs", "control_agents", "background_agents", "contrast_reversed", "is_combination", "dose_nM",
         "inhibitor_primary_targets", "inhibitor_primary_targets_in_universe", "inhibitor_profiled_targets", "target_confidence",
         "genetic_target_kinases", "stimulus", "stimulus_class", "stimulus_expected_up_kinases",
         "cell_type", "cell_line_norm", "sample_type", "time_raw", "time_min", "time_source", "time_note",
         "time_series_id", "time_series_n_timepoints", "time_rank", "treatment_key", "dose_series_id",
         "labeling", "spectrometer", "n_sites", "n_psp_ks_sites", "qc_flag", "duplicate_of"]
A = A[order]
A.to_csv(RES / "conditions_annotated.tsv", sep="\t", index=False)
cnt = A.groupby(["perturbation_class"]).agg(n_conditions=("conditionID", "size"), n_studies=("study_pmid", "nunique")).reset_index()
cnt.to_csv(RES / "condition_class_counts.tsv", sep="\t", index=False)
print(cnt.to_string())
print(A["perturbation_subclass"].value_counts().head(40).to_string())
print("time series:", A["time_series_id"].ne("").sum(), "conditions in", A.loc[A.time_series_id.ne(""), "time_series_id"].nunique(), "series")
print("dose series:", A["dose_series_id"].ne("").sum(), "conditions in", A.loc[A.dose_series_id.ne(""), "dose_series_id"].nunique())
print("time_min parsed:", A["time_min"].notna().sum())
print("target conf:", A["target_confidence"].value_counts().to_dict())
print("qc:", A["qc_flag"].value_counts().to_dict())
# agreement with the kinaseActivityMain parser
dd = dose[dose.side == "condition"].groupby("conditionID")["agent"].agg(lambda x: set(x))
mine = A.set_index("conditionID")["agents"].map(lambda s: set(filter(None, s.split(";"))))
cmp = pd.concat([dd.rename("kam"), mine.rename("mine")], axis=1)
cmp["kam"] = cmp["kam"].map(lambda x: x if isinstance(x, set) else set())
dis = cmp[cmp["kam"] != cmp["mine"]]
dis.to_csv(RES / "_agent_parse_disagreements_vs_kam.tsv", sep="\t")
print("agent-parse disagreements vs kinaseActivityMain condition_doses:", len(dis))
