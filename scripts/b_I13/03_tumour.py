"""Task 2/3: phosphatase lesions in CPTAC tumours as natural experiments.
Stage 1 (per cohort): site z (protein-adjusted residual, z-scored, NaN->0, sites >=50% observed) ~ exposure + purity + E2F + G2M + CNV burden + immune/stromal (FWL, vectorised).
Stage 2: per exposure, site partial-beta ~ phosphatase-preference feature (+ site nuisance). Null = same pipeline on frequency-matched random genes (other chromosomes) and same-arm genes.
"""
import sys, re
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/scripts/infra_cptac"); import cptac_load as cl
NPC = int(sys.argv[1]) if len(sys.argv) > 1 else 0
TAG = f"_pc{NPC}" if NPC else ""
R = "/workspace/results/b_I13"; rng = np.random.default_rng(13)
F = pd.read_parquet(f"{R}/cache/site_features.parquet").set_index("site")
W = np.array([list(w) for w in F.window])
F["akt_motif"] = ((W[:, 2] == "R") & (W[:, 4] == "R")).astype(int)   # R-x-R-x-x-S/T
F["SQ_acid"] = ((F.Q1 == 1) & (F.acid_m7p7 >= 2)).astype(int)
HG = pd.read_csv("<HACKATHON_ROOT>/kinaseSubMain/external/hgnc/hgnc_complete_set.txt", sep="\t", usecols=["symbol", "location"], low_memory=False).dropna()
HG["arm"] = HG.location.str.extract(r"^(\w+?[pq])")[0]; ARM = HG.set_index("symbol").arm
UNI = {"PPP2R2A": "P63151", "PPP2R1A": "P30153", "PPM1D": "O15297", "PTPRD": "P23468", "PTPN12": "Q05209", "PTPN11": "Q06124", "PTPN1": "P18031", "PTEN": "P60484", "PPP6C": "O00743"}
HOT = {"PPP2R1A": r"p\.(P179|R182|R183|S256|W257|R258)"}

def covars(c, t):
    cases = [s[:-2] for s in t]
    cli = cl.table(c, "clinical").drop_duplicates("case_id").set_index("case_id")
    ann = cl.table(c, "annotations").drop_duplicates("case_id").set_index("case_id")
    cnv = cl.matrix(c, "cnv_log2")
    burden = cnv.abs().mean(0)
    Z = pd.DataFrame(index=cases)
    Z["purity"] = cli.WES_purity.reindex(cases).values
    for k in ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore"]:
        Z[k] = ann[k].reindex(cases).values if k in ann else np.nan
    Z["cnv_burden"] = burden.reindex(cases).values
    Z = Z.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0))
    Z = (Z - Z.mean()) / Z.std().replace(0, 1); Z = Z.fillna(0)
    Z.index = t; return Z

def ntype(e): return "null_prot" if "_prot" in e else ("null_perm" if "hotspot" in e else "null_other")

def resid(M, Q):  # residualise columns of M on orthonormal Q
    return M - Q @ (Q.T @ M)

stage1, stage2, firststage, specific, kinact = [], [], [], [], []
SPEC = {"TP53_S15": "P04637_S15", "CHEK1_S345": "O14757_S345", "CHEK1_S317": "O14757_S317", "H2AX_S140": "P16104_S140", "ATM_S1981": "Q13315_S1981",
        "AKT1_S473": "P31749_S473", "AKT1_T308": "P31749_T308", "AKT1S1_T246": "Q96B36_T246", "GSK3B_S9": "P49841_S9", "STAT3_Y705": "P40763_Y705",
        "MDM2_S166": "Q00987_S166", "KAP1_S824": "Q13263_S824", "SMC1A_S957": "Q14683_S957", "53BP1_S1778": "Q12888_S1778"}
COMMON = None
for c in cl.COHORTS:
    p = cl.matrix(c, "protein"); t_ = [x for x in cl.tumour_ids(c) if x in p.columns]
    ok = set(p.index[p[t_].notna().mean(1) >= 0.9]); COMMON = ok if COMMON is None else COMMON & ok
NULLPROT = sorted(rng.choice(sorted(COMMON), 1000, replace=False)); print("common proteins", len(COMMON))
KZ = {v: cl.activity(v, "corev3_z") for v in ["tc", "tc_padj"]}
for c in cl.COHORTS:
    t = cl.tumour_ids(c); cases = [s[:-2] for s in t]
    ph = cl.matrix(c, "input_tc_padj"); ph = ph[[s for s in t if s in ph.columns]]; t = list(ph.columns); cases = [s[:-2] for s in t]
    ph = ph[ph.notna().mean(1) >= 0.5]
    Y = ph.sub(ph.mean(1), axis=0).div(ph.std(1), axis=0).fillna(0).T.values  # samples x sites
    sites = ph.index
    Zc = covars(c, t)
    if NPC > 0:   # latent proteome factors (top PCs of protein matrix, >=90% observed proteins)
        pp = cl.matrix(c, "protein").reindex(columns=t).astype(float); pp = pp[pp.notna().mean(1) >= 0.9]
        pp = pp.sub(pp.mean(1), axis=0).T.fillna(0).values
        U, sv, _ = np.linalg.svd(pp - pp.mean(0), full_matrices=False)
        for i in range(NPC): Zc[f"protPC{i+1}"] = U[:, i]
    Q, _ = np.linalg.qr(np.column_stack([np.ones(len(t)), Zc.values]))
    Yr = resid(Y, Q)
    cnv = cl.matrix(c, "cnv_log2"); gis = cl.matrix(c, "cnv_gistic"); prot = cl.matrix(c, "protein"); mut = cl.table(c, "mutations")
    cnv = cnv[~cnv.index.duplicated()].reindex(columns=cases); gis = gis[~gis.index.duplicated()].reindex(columns=cases).astype(float); prot = prot[~prot.index.duplicated()]
    ok = cnv.notna().mean(1) > 0.9; cnv = cnv[ok].apply(lambda s: s.fillna(s.median()), axis=1)
    # exposures (oriented so positive = less phosphatase / for PPM1D: more phosphatase)
    EX = {}
    for g in ["PPP2R2A", "PTPRD", "PTPN12", "PTPN11", "PTPN1", "PTEN", "PPP6C"]:
        if g in cnv.index: EX[f"{g}_cnvloss"] = -cnv.loc[g].values
        if UNI[g] in prot.index:
            pv = prot.loc[UNI[g]].reindex(t)
            if pv.notna().mean() > 0.8: EX[f"{g}_protlow"] = -pv.fillna(pv.median()).values
        if g in cnv.index and UNI[g] in prot.index:
            pv = prot.loc[UNI[g]].reindex(t).values.astype(float); cv = cnv.loc[g].values; m = ~np.isnan(pv)
            if m.sum() > 20: firststage.append(dict(cohort=c, gene=g, rho=pd.Series(cv[m]).corr(pd.Series(pv[m]), method="spearman"), n=int(m.sum()),
                                                  prot_shallowdel_minus_neutral=np.nanmean(pv[(gis.loc[g].values <= -1) & m]) - np.nanmean(pv[(gis.loc[g].values == 0) & m]) if g in gis.index else np.nan,
                                                  n_del=int((gis.loc[g].values <= -1).sum()) if g in gis.index else 0))
    pv = prot.loc[UNI["PPP2R2A"]].reindex(t).astype(float); pv = pv.fillna(pv.median()).values; cv = cnv.loc["PPP2R2A"].values
    A = np.column_stack([np.ones(len(cv)), cv]); EX["PPP2R2A_protlow_residCNV"] = -(pv - A @ np.linalg.lstsq(A, pv, rcond=None)[0])
    if "PPM1D" in cnv.index:
        EX["PPM1D_cnvgain"] = cnv.loc["PPM1D"].values
        pv = prot.loc[UNI["PPM1D"]].reindex(t) if UNI["PPM1D"] in prot.index else None
        if pv is not None and pv.notna().mean() > 0.5: EX["PPM1D_prothigh"] = pv.fillna(pv.median()).values
        if pv is not None: 
            m = pv.notna().values
            firststage.append(dict(cohort=c, gene="PPM1D", rho=pd.Series(cnv.loc["PPM1D"].values[m]).corr(pd.Series(pv.values[m]), method="spearman"), n=int(m.sum())))
    EXTRA = {"PPP2CA": "P67775", "PPP2R1A": "P30153", "PPP2R5A": "Q15172", "PPP2R5C": "Q13362", "PPP2R5D": "Q14738", "PPP2R5E": "Q16537", "PPP2R2D": "Q66LE6",
             "PPP1CA": "P62136", "PPP1CB": "P62140", "PPP1CC": "P36873", "MASTL": "Q96GX5", "PTPN2": "P17706", "PTPN6": "P29350"}
    for g, u in EXTRA.items():
        if u in prot.index:
            pv = prot.loc[u].reindex(t).astype(float)
            if pv.notna().mean() > 0.8: EX[f"{g}_protlow"] = -pv.fillna(pv.median()).values
    for g, u in {"ARPP19": "P56211", "ENSA": "O43768"}.items():   # B55 inhibitors: high = less B55 activity
        if u in prot.index:
            pv = prot.loc[u].reindex(t).astype(float)
            if pv.notna().mean() > 0.8: EX[f"{g}_prothigh"] = pv.fillna(pv.median()).values
    hot = mut[(mut.gene == "PPP2R1A") & mut.protein_change.str.contains(HOT["PPP2R1A"], na=False)].case_id.unique()
    hv = np.isin(cases, hot).astype(float)
    if hv.sum() >= 3: EX["PPP2R1A_hotspot"] = hv
    # TP53 mutation as extra covariate for PPM1D exposures is handled by including it in exposure residualisation
    tp53 = np.isin(cases, mut[mut.gene == "TP53"].case_id.unique()).astype(float)
    names = list(EX); E = np.column_stack([EX[k] for k in names])
    # null genes: (a) other-chromosome genes, (b) same-arm (8p) genes
    gl = [g for g in cnv.index if g in ARM.index]
    other = [g for g in gl if not str(ARM[g]).startswith("8")]
    arm8p = [g for g in gl if ARM[g] == "8p" and g != "PPP2R2A"]
    nul_o = list(rng.choice(other, 300, replace=False)); nul_8p = list(rng.choice(arm8p, min(60, len(arm8p)), replace=False))
    En = -cnv.loc[nul_o + nul_8p].values.T
    pg = prot.reindex(columns=t).astype(float); pg = pg[pg.notna().mean(1) >= 0.9]
    nul_p = [u for u in NULLPROT if u in pg.index]
    Ep = -pg.loc[nul_p].T.apply(lambda s: s.fillna(s.median())).values
    nperm = []
    if "PPP2R1A_hotspot" in EX:
        k = int(EX["PPP2R1A_hotspot"].sum())
        for i in range(300):
            v = np.zeros(len(t)); v[rng.choice(len(t), k, replace=False)] = 1; nperm.append(v)
    En = np.column_stack([En, Ep] + ([np.column_stack(nperm)] if nperm else []))
    Eall = np.column_stack([E, En]); Q2, _ = np.linalg.qr(np.column_stack([Q, tp53]))
    Er = resid(Eall, Q2); Yr2 = resid(Y, Q2)
    B = (Er.T @ Yr2) / (Er ** 2).sum(0)[:, None]               # exposure x site partial slope (site z per unit exposure)
    sdE = Er.std(0); B = B * sdE[:, None]                       # per SD of exposure
    lab = names + [f"null_other:{g}" for g in nul_o] + [f"null_8p:{g}" for g in nul_8p] + [f"null_prot:{g}" for g in nul_p] + [f"null_perm:{i}" for i in range(len(nperm))]
    Bdf = pd.DataFrame(B, index=lab, columns=sites)
    Bdf.loc[names].T.astype("float32").to_parquet(f"{R}/cache/stage1_{c}{TAG}.parquet")
    for k, sid in SPEC.items():
        if sid in Bdf.columns:
            for e in names:
                nulls = Bdf.loc[[l for l in lab if l.startswith(ntype(e))], sid]
                specific.append(dict(cohort=c, exposure=e, site=k, beta=Bdf.loc[e, sid], null_mean=nulls.mean(), null_sd=nulls.std(), z_emp=(Bdf.loc[e, sid] - nulls.mean()) / nulls.std()))
    # stage 2
    Fs = F.reindex(sites)
    raw = cl.matrix(c, "phospho").reindex(index=sites, columns=t)
    nuis = pd.DataFrame({"sd": np.log(raw.std(1).values + 1e-3), "mean": raw.mean(1).values, "obs": ph.notna().mean(1).values}, index=sites)
    nuis = nuis.fillna(nuis.median())
    TESTS = {"B55_lit|ST": ("B55_lit", Fs.residue.isin(["S", "T"])), "isT|ST": ("isT", Fs.residue.isin(["S", "T"])),
             "P1|ST": ("P1", Fs.residue.isin(["S", "T"])), "bipartite_basic|ST": ("bipartite_basic", Fs.residue.isin(["S", "T"])),
             "basic_m3p3|ST": ("basic_m3p3", Fs.residue.isin(["S", "T"])), "acid_m3p3|ST": ("acid_m3p3", Fs.residue.isin(["S", "T"])),
             "B56_lit|ST": ("B56_lit", Fs.residue.isin(["S", "T"])), "R_m3|ST": ("R_m3", Fs.residue.isin(["S", "T"])),
             "isY|all": ("isY", Fs.residue.notna()), "Q1|ST": ("Q1", Fs.residue.isin(["S", "T"])), "SQ_acid|ST": ("SQ_acid", Fs.residue.isin(["S", "T"])),
             "akt_motif|ST": ("akt_motif", Fs.residue.isin(["S", "T"]))}
    for tn, (feat, mask) in TESTS.items():
        m = mask.fillna(False).values & Fs[feat].notna().values
        x = Fs[feat].values[m].astype(float); x = (x - x.mean()) / x.std()
        Xn = np.column_stack([np.ones(m.sum()), nuis.values[m]]); Qn, _ = np.linalg.qr(Xn)
        xr = x - Qn @ (Qn.T @ x)
        Bm = B[:, m]; Bmr = Bm - (Bm @ Qn) @ Qn.T
        coef = (Bmr @ xr) / (xr @ xr)                                   # per exposure
        res = Bmr - np.outer(coef, xr); se = np.sqrt((res ** 2).sum(1) / (m.sum() - Xn.shape[1] - 1) / (xr @ xr))
        for i, l in enumerate(lab):
            stage2.append(dict(cohort=c, test=tn, exposure=l, coef=coef[i], se_naive=se[i], n_sites=int(m.sum()), n_tum=len(t),
                               n_exposed=int((EX[l] > 0).sum()) if l == "PPP2R1A_hotspot" else np.nan))
    # Task 3: kinase-activity bias
    for v, Kz in KZ.items():
        Kc = Kz.reindex(columns=t).T.fillna(0).values
        Kr = resid(Kc, Q2); KB = (Er.T @ Kr) / (Er ** 2).sum(0)[:, None] * sdE[:, None] / Kr.std(0).clip(1e-6)[None, :]
        kd = pd.DataFrame(KB, index=lab, columns=Kz.index)
        for kin in ["CDK1", "CDK2", "CDK5", "MAPK1", "MAPK3", "GSK3B", "AKT1", "PRKACA", "AURKA", "AURKB", "PLK1", "CSNK2A1", "ATM", "ATR", "CHEK1", "PRKDC", "SRC", "EGFR", "MTOR", "HIPK2", "DYRK1A", "CDK9", "MAPK14", "RPS6KB1"]:
            if kin not in kd.columns: continue
            for e in names:
                nn = kd.loc[[l for l in lab if l.startswith(ntype(e))], kin]
                n8 = kd.loc[[l for l in lab if l.startswith("null_8p")], kin]
                kinact.append(dict(cohort=c, variant=v, exposure=e, kinase=kin, beta=kd.loc[e, kin], null_mean=nn.mean(), null_sd=nn.std(), null8p_mean=n8.mean(), null8p_sd=n8.std(), n_null8p=len(n8)))
    print(c, len(t), len(sites), names, flush=True)
pd.DataFrame(stage2).to_parquet(f"{R}/cache/stage2_all{TAG}.parquet")
pd.DataFrame(firststage).to_csv(f"{R}/tables/t2_first_stage_cnv_protein{TAG}.tsv", sep="\t", index=False)
pd.DataFrame(specific).to_csv(f"{R}/tables/t2_specific_sites{TAG}.tsv", sep="\t", index=False)
pd.DataFrame(kinact).to_csv(f"{R}/tables/t3_kinase_activity_vs_lesion{TAG}.tsv", sep="\t", index=False)
