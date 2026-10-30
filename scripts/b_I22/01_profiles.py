"""Per-experiment primary/reverse kinase-motif profiles (inh + stim), two definitions, raw and study x cell centred.
lor90: I05 definition (responders |x|>1 in direction, >=20; bg |x|<0.3; log OR of Cantley pct>=90).
dpct : mean Cantley percentile among responders (|x|>0.58, >=30) minus mean among bg (|x|<0.3) — less noisy.
Also stores measured-site masks (for coverage Jaccard) and per-experiment metadata."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import load_exp, mask_binders, CA as ACA
OUT = '/workspace/cache/b_I22/'
MP = pd.read_parquet('/workspace/results/b_I05/cache/motif_pct.parquet')
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'residue'], low_memory=False).set_index('site')
KST = pd.read_csv('<HACKATHON_ROOT>/kinaseSubMain/external/kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases/kinase_data/kinome_information.tsv', sep='\t')
kst = [k for k in KST[KST.TYPE == 'ser_thr'].GENE_NAME.unique() if k in MP.columns]

def get(kind):
    if kind == 'inh':
        X, E = load_exp('inh'); A = pd.DataFrame(mask_binders(X, E), index=X.index, columns=X.columns)
        E = E.rename(columns={E.columns[0]: 'exp'}); E = E[['exp', 'study', 'cell', 'cls', 'agent', 'time_min', 'dose_nM']].copy()
    else:
        M = pd.read_csv(ACA + 'stim_meta.tsv', sep='\t', dtype={'study_pmid': str})
        drop = (M.control_raw == 'LPA') | (M.study_pmid == '29020638') | M.condition_raw.str.contains('p38i', na=False) | M.agent.isin(['cytarabine', 'methotrexate'])
        M = M[~drop]; Xc = pd.read_parquet(ACA + 'stim_cond.parquet')[M.conditionID]
        M['exp'] = M.study_pmid + '|' + M.cell_line_norm.astype(str) + '|' + M.cls
        A = Xc.T.groupby(M.exp.values).mean().T
        E = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), cls=('cls', 'first'), agent=('agent', 'first'),
                                 time_min=('time_min', 'median')).reindex(A.columns).reset_index()
        E['dose_nM'] = np.nan
    A = A[S.residue.reindex(A.index).isin(['S', 'T']).values]
    E['kind'] = kind; E['cell'] = E.cell.astype(str)
    return A, E.reset_index(drop=True)

def centre(A, E):
    A = A.copy(); g = E.study + '|' + E.cell
    for k, d in E.groupby(g):
        if len(d) >= 5:
            cols = A.columns[d.index.values]
            A[cols] = A[cols].sub(A[cols].median(1), axis=0)
    return A

profs = {}; metas = []; masks = []
for kind in ['inh', 'stim']:
    A, E = get(kind)
    masks.append(A.notna())
    for cen in ['raw', 'cen']:
        B = centre(A, E) if cen == 'cen' else A
        P = MP.reindex(B.index)[kst].values.astype(np.float32)
        M90 = P >= 90
        for j in range(B.shape[1]):
            x = B.iloc[:, j].values; fin = ~np.isnan(x)
            if fin.sum() < 500: continue
            xx = np.nan_to_num(x); bg = fin & (np.abs(xx) < 0.3)
            prim = -1 if kind == 'inh' else 1
            for lab, sg in [('primary', prim), ('reverse', -prim)]:
                for defn, thr, nmin in [('lor90', 1, 20), ('dpct', 0.58, 30)]:
                    rsp = fin & (sg * xx > thr)
                    if rsp.sum() < nmin: continue
                    if defn == 'lor90':
                        a = M90[rsp].sum(0) + .5; b = rsp.sum() - M90[rsp].sum(0) + .5; c = M90[bg].sum(0) + .5; d = bg.sum() - M90[bg].sum(0) + .5
                        v = np.log(a * d / (b * c))
                    else:
                        v = np.nanmean(P[rsp], 0) - np.nanmean(P[bg], 0)
                    profs.setdefault((defn, cen), []).append(v)
                    r = E.iloc[j].to_dict(); r.update(dir=lab, defn=defn, cen=cen, n_resp=int(rsp.sum()), n_meas=int(fin.sum()))
                    metas.append(r)
Mt = pd.DataFrame(metas)
allP = np.vstack([np.vstack(profs[k]) for k in profs]); 
# order metas the same as stacking: rebuild
order = []
for k in profs: order.append(Mt[(Mt.defn == k[0]) & (Mt.cen == k[1])])
Mt = pd.concat(order).reset_index(drop=True)
pd.DataFrame(allP, columns=kst).to_parquet(OUT + 'profiles.parquet'); Mt.to_csv(OUT + 'profiles_meta.tsv', sep='\t', index=False)
masks[0].columns = ['inh::' + c for c in masks[0].columns]; masks[1].columns = ['stim::' + c for c in masks[1].columns]
Mask = pd.concat(masks, axis=1)
Mask.fillna(False).astype(np.uint8).to_parquet(OUT + 'measured_mask.parquet')
print(Mt.groupby(['defn', 'cen', 'kind', 'dir']).size())
