"""I26 step 6: summarise activity benchmark (paired vs base; study-clustered)."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
sys.path.insert(0, os.path.dirname(__file__)); from lib26 import *

R = pd.read_parquet(f"{C26}/activity_benchmark_long.parquet")
R["vclass"] = R.variant.str.replace(r"_\d+$", "", regex=True)
R["sch"] = np.where(R.scheme == "LOSO", "LOSO", "half-split x20")
# average random reps per condition
A = R.groupby(["sch", "scheme", "condition", "study", "target_confidence", "vclass"]).agg(auroc=("auroc", "mean"), tz=("mean_target_z", "mean")).reset_index()
W = A.pivot_table(index=["sch", "scheme", "condition", "study", "target_confidence"], columns="vclass", values=["auroc", "tz"]).reset_index()
W.columns = ["_".join([c for c in col if c]) for col in W.columns]
W.to_csv(f"{R26}/tables/activity_benchmark_conditions.tsv", sep="\t", index=False, float_format="%.4g")
rows = []
for sch, d0 in W.groupby("sch"):
    for conf in ["high+medium", "high"]:
        d1 = d0 if conf == "high+medium" else d0[d0.target_confidence == "high"]
        for sub in ["all", "affected"]:
            d = d1 if sub == "all" else d1[(d1.auroc_prune - d1.auroc_base).abs() > 1e-9]
            d = d.dropna(subset=["auroc_base", "auroc_prune"])
            for v in ["prune", "prune+amb", "prune+add", "rand", "rand_any"]:
                dd = d.dropna(subset=[f"auroc_{v}"])
                delta = dd[f"auroc_{v}"] - dd.auroc_base
                dtz = dd[f"tz_{v}"] - dd.tz_base
                # study-clustered: per (scheme, study) mean delta
                st = delta.groupby([dd.scheme, dd.study]).mean()
                try:
                    p = wilcoxon(st[st != 0]).pvalue if (st != 0).sum() >= 5 else np.nan
                except ValueError:
                    p = np.nan
                rows.append(dict(scheme=sch, targets=conf, subset=sub, variant=v, n_cond_evals=len(dd), n_studies=dd.study.nunique(),
                                 median_auroc_base=dd.auroc_base.median(), median_auroc_var=dd[f"auroc_{v}"].median(),
                                 mean_auroc_base=dd.auroc_base.mean(), mean_auroc_var=dd[f"auroc_{v}"].mean(),
                                 mean_delta=delta.mean(), frac_improved=(delta > 1e-9).mean(), frac_worse=(delta < -1e-9).mean(),
                                 study_mean_delta=st.mean(), study_frac_improved=(st > 0).mean(), wilcoxon_p_studies=p,
                                 mean_delta_target_z=dtz.mean()))
S = pd.DataFrame(rows)
S.to_csv(f"{R26}/tables/activity_benchmark_summary.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(S.round(4).to_string())
# per study LOSO
d = W[W.sch == "LOSO"]
print((d.assign(dp=d["auroc_prune"] - d.auroc_base, dr=d["auroc_rand"] - d.auroc_base, da=d["auroc_prune+add"] - d.auroc_base)
       .groupby("study")[["auroc_base", "dp", "dr", "da"]].agg(["mean", "size"]).round(3).to_string()))
