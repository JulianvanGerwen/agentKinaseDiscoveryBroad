"""Kd-independent amplification index: potency of sites ON the target proteins (target cluster T, 'self';
for ABL1 also BCR = BCR-ABL fusion) vs (a) a pre-specified panel of canonical downstream readouts and
(b) all other down-regulated sites in the same experiment.
amp_index = median pEC50(readouts) - median pEC50(self down-sites).  >0 = downstream responds at lower dose
than target autophosphorylation (signal amplification / receptor reserve).
Null: within-experiment label permutation of self vs readout; sign test across drugs.
Outputs: results/b_I17/amplification_by_experiment.tsv, amplification_sites.tsv
"""
import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
RES = "/workspace/results/b_I17"
A = pd.read_csv(f"{RES}/shift_sites.tsv", sep="\t")
A["residue"] = A.site.str.extract(r"_([STY])\d+$")[0]
READOUT = {"MAPK1": ["T185", "Y187"], "MAPK3": ["T202", "Y204"], "RPS6KA1": ["S380", "T573", "S221", "T359", "S363"],
           "RPS6KA3": ["S227", "S369", "T577", "Y580", "S386"], "AKT1S1": ["T246"], "GSK3B": ["S9"], "GSK3A": ["S21"],
           "STAT5A": ["Y694"], "STAT5B": ["Y699"], "CRKL": ["Y207"], "RPS6": ["S235", "S236", "S240", "S244"],
           "EIF4EBP1": ["T37", "T46", "S65", "T70"], "MAPKAPK2": ["T334"], "PTK2": ["Y397", "Y576"], "PXN": ["Y118", "Y31"]}
A["pos"] = A.site.str.split("_").str[1]
A["readout"] = [g in READOUT and p in READOUT[g] for g, p in zip(A.gene, A.pos)]
E = pd.read_csv("/workspace/cache/b_I17/experiments.tsv", sep="\t").set_index("exp_id")
prim = {ex: set(str(E.loc[ex, "targets"]).split(";")) for ex in E.index}
A["self2"] = A.self | ((A.target == "ABL1") & (A.gene == "BCR")) | np.array([g in prim.get(ex, set()) for g, ex in zip(A.gene, A.exp_id)])
dn = A[A.direction == "down"].copy()
rng = np.random.default_rng(11)
rows = []
for ex, d in dn.groupby("exp_id"):
    s = d[d.self2]; r = d[d.readout & ~d.self2]
    pct_self = np.mean([(d.pEC50 > v).mean() for v in s.pEC50]) if len(s) else np.nan  # fraction of down sites more potent than a self site
    row = dict(exp_id=ex, agent=d.agent.iloc[0], cell_line=d.cell_line.iloc[0], target=d.target.iloc[0], pKd_target=d.pKd_target.iloc[0],
               n_down=len(d), n_self=len(s), n_readout=len(r), med_self=s.pEC50.median(), med_readout=r.pEC50.median(),
               med_all_down=d.pEC50.median(), frac_down_more_potent_than_self=pct_self,
               hill_self=s.hill.median(), hill_readout=r.hill.median())
    if len(s) >= 2 and len(r) >= 2:
        obs = r.pEC50.median() - s.pEC50.median()
        both = pd.concat([s.pEC50, r.pEC50]).values; ns = len(s)
        null = []
        for _ in range(2000):
            p = rng.permutation(both); null.append(np.median(p[ns:]) - np.median(p[:ns]))
        row.update(amp_index=obs, p_perm=(np.sum(np.array(null) >= obs) + 1) / 2001)
    rows.append(row)
R = pd.DataFrame(rows)
R.to_csv(f"{RES}/amplification_by_experiment.tsv", sep="\t", index=False)
dn[dn.self2 | dn.readout].to_csv(f"{RES}/amplification_sites.tsv", sep="\t", index=False)
print(R.round(2).to_string())
v = R.dropna(subset=["amp_index"])
dl = v.groupby("agent").amp_index.median()
print("experiments with amp_index>0:", (v.amp_index > 0).sum(), "/", len(v), "| drugs:", (dl > 0).sum(), "/", len(dl),
      "binom p (drug level)", binomtest(int((dl > 0).sum()), len(dl), 0.5, alternative="greater").pvalue,
      "| median amp (drug level)", dl.median(), "| wilcoxon drug-level p", wilcoxon(dl, alternative="greater").pvalue if len(dl) > 5 else None)
w = R.dropna(subset=["frac_down_more_potent_than_self"])
print("median fraction of down-sites more potent than self sites:", w.frac_down_more_potent_than_self.median(), "n exp", len(w))
