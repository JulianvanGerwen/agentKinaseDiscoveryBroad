"""b_Q8 step 8: diagnostics. (a) How much of each footprint score's variance lies in the top-3 global PCs, in NAT vs tumour
(explains why PC adjustment changes the sign of MAPK/RTK Delta). (b) Gene-annotated per-cohort detail of robust rewired site edges."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from common import *
rows = []
FPs = ['KI_MEK', 'EGF', 'GF_all', 'KI_MET_ALK', 'KI_EGFR', 'GPCR_ligand', 'DNA_damage', 'metformin', 'TNF']
for coh in MAINC:
    ph, pr, T, N, paired = load_raw(coh)
    cs = sorted(c for c in paired if f'{c}_N' in N and f'{c}_T' in T)
    for grp, smp in [('T', [c + '_T' for c in cs]), ('N', [c + '_N' for c in cs])]:
        Z = process(ph, pr, smp, padj=True, k=0, minobs=max(8, int(0.4 * len(smp))))
        F, ve = global_pcs(Z, 3)
        for g in FPs:
            s = score(Z, fp_set(g)); m = s.notna().values
            X = np.column_stack([np.ones(m.sum()), F[m]]); b = np.linalg.lstsq(X, s.values[m], rcond=None)[0]
            r2 = 1 - ((s.values[m] - X @ b) ** 2).sum() / ((s.values[m] - s.values[m].mean()) ** 2).sum()
            rows.append(dict(cohort=coh, grp=grp, fp=g, R2_top3PC=r2, PCvar3=ve[:3].sum()))
D = pd.DataFrame(rows); D.to_csv(TQ + 'T10_fp_variance_in_global_PCs.tsv', sep='\t', index=False)
print(D.pivot_table(index='fp', columns='grp', values='R2_top3PC', aggfunc='median').round(2))
print(D.groupby('grp').PCvar3.median().round(3))
# (b)
R = pd.read_csv(TQ + 'T8_edge_rewiring.tsv', sep='\t', index_col=0)
D0 = pd.read_parquet(CQ + 'edge_stats_padj_pc3.parquet')
S = site_annot()
rw = R[R.rewired].copy()
rw['site'] = [r.split('>')[1] if '>' in r else '' for r in rw.index]
rw['gene_site'] = [f"{S.gene.get(s, '?')}_{s.split('_')[1]}" if s else '' for s in rw.site]
det = D0[D0.rel.isin(rw.index) & D0.cohort.isin(MAINC)].pivot(index='rel', columns='cohort', values=['r_N', 'r_T'])
det.columns = [f'{a}_{b}' for a, b in det.columns]
rw = rw.join(det.round(2)); rw.to_csv(TQ + 'T8b_rewired_detail.tsv', sep='\t')
print(rw.sort_values('meta_z')[['gene_site', 'rN_med', 'rT_med', 'meta_z', 'meta_q', 'n_up', 'n_down', 'meta_z_padj', 'meta_z_comp_pc3']].round(2).to_string())
