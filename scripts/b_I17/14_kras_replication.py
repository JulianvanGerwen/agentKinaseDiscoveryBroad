"""Replication on an independent Kuster-lab dataset (Zenodo 16096119, 'Illuminating oncogenic KRAS signaling...',
CurveCurator phospho fits, 2 h, 3 cell lines). decryptM 2.0 (Zenodo 17533475) is access-restricted (HTTP 403).
Kinase inhibitors present: trametinib (allosteric MEK, 3 lines), temuterkib (ATP-site ERK1/2, 2 lines; no Klaeger/Davis Kd).
Readouts: ERK1/2 activation-loop peptides (direct MEK substrates), and all significantly down-regulated phosphopeptides.
Shift vs Klaeger trametinib pKd (MAP2K2 9 nM / MAP2K1 62 nM). Output: results/b_I17/kras_replication.tsv"""
import zipfile, re, numpy as np, pandas as pd
Z = zipfile.ZipFile("/workspace/cache/b_I17/kras_zenodo/CurveCurator_Phospho.zip")
files = [n for n in Z.namelist() if re.search(r"/curves_[^/]*\.txt$", n)]
ERK = {"MAPK1": "VADPDHDHTGFLTEYVATR", "MAPK3": "IADPEHDHTGFLTEYVATR"}
rows, rd = [], []
for f in files:
    exp = f.split("/")[1]
    if not any(k in exp for k in ["Trametinib", "Temuterkib", "RMC4630", "BI3406", "Sotorasib_2h_R1", "Adagrasib", "MRTX1133"]): continue
    d = pd.read_csv(Z.open(f), sep="\t", low_memory=False)
    seq = d["Modified sequence"].astype(str).str.replace(r"\(.*?\)", "", regex=True).str.replace("_", "")
    sig_dn = (d["Curve Regulation"] == "down")
    for g, pep in ERK.items():
        m = seq.str.contains(pep) & sig_dn
        for _, r in d[m].iterrows():
            rd.append(dict(exp=exp, gene=g, peptide=r["Modified sequence"], pEC50=r.pEC50, slope=r["Curve Slope"], fc=r["Curve Fold Change"]))
    rows.append(dict(exp=exp, n_peptides=len(d), n_down=int(sig_dn.sum()), n_up=int((d["Curve Regulation"] == "up").sum()),
                     median_pEC50_down=d.loc[sig_dn, "pEC50"].median(), median_slope_down=d.loc[sig_dn, "Curve Slope"].median(),
                     ERK_pEC50=np.median([x["pEC50"] for x in rd if x["exp"] == exp]) if any(x["exp"] == exp for x in rd) else np.nan))
R = pd.DataFrame(rows)
pkd_mek2, pkd_mek1 = 9 - np.log10(9), 9 - np.log10(62)
R["shift_ERK_vs_MEK2"] = np.where(R.exp.str.contains("Trametinib"), R.ERK_pEC50 - pkd_mek2, np.nan)
R["shift_all_vs_MEK2"] = np.where(R.exp.str.contains("Trametinib"), R.median_pEC50_down - pkd_mek2, np.nan)
R.to_csv("/workspace/results/b_I17/kras_replication.tsv", sep="\t", index=False)
pd.DataFrame(rd).to_csv("/workspace/results/b_I17/kras_replication_erk_peptides.tsv", sep="\t", index=False)
print(R.round(2).to_string()); print(pd.DataFrame(rd).round(2).to_string())
