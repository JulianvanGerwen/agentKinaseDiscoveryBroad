"""Cheng-Prusoff consistency with ATP-free Davis 2011 Kd: implied cellular [ATP] = Km*(10^-shift - 1) per ATP-site drug.
Compares dispersion of log implied [ATP] (Km-aware model) vs dispersion of the raw shift (constant-offset model);
null = permute Km over drugs (all 5! orders) ; LODO of median implied ATP. Output: results/b_I17/atp_implied.tsv"""
import numpy as np, pandas as pd, itertools
RES = "/workspace/results/b_I17"
D = pd.read_csv(f"{RES}/atp_drug_table.tsv", sep="\t")
d = D[D.atp_site & D.Km_ATP_uM.notna() & D.shift_tsub_davis.notna()].copy()
d["implied_ATP_mM"] = d.Km_ATP_uM * (10 ** (-d.shift_tsub_davis) - 1) / 1000
d["implied_ATP_mM_kinobeads"] = d.Km_ATP_uM * (10 ** (-d.shift_tsub) - 1) / 1000
obs_sd = np.log10(d.implied_ATP_mM).std(ddof=1); off_sd = d.shift_tsub_davis.std(ddof=1)
null = []
for perm in itertools.permutations(d.Km_ATP_uM.values):
    null.append(np.log10(np.array(perm) * (10 ** (-d.shift_tsub_davis.values) - 1)).std(ddof=1))
null = np.array(null)
lodo = [np.median(d.implied_ATP_mM.drop(i)) for i in d.index]
d.to_csv(f"{RES}/atp_implied.tsv", sep="\t", index=False)
print(d[["agent", "target", "Km_ATP_uM", "davis_kd_nM", "shift_tsub_davis", "implied_ATP_mM", "shift_tsub", "implied_ATP_mM_kinobeads"]].round(3).to_string())
print(f"median implied ATP {d.implied_ATP_mM.median():.2f} mM (range {d.implied_ATP_mM.min():.2f}-{d.implied_ATP_mM.max():.2f}); LODO median {min(lodo):.2f}-{max(lodo):.2f}")
print(f"SD log10 implied ATP {obs_sd:.3f} vs SD shift (constant offset) {off_sd:.3f}; Km-permutation p = {(null <= obs_sd + 1e-12).mean():.3f} (n perms {len(null)})")
print("non-covalent only SD:", np.log10(d[~d.covalent].implied_ATP_mM).std(ddof=1))
