"""Task 3b: which kinase motifs are enriched among reverse-direction responders, and are the profiles consistent across
contexts (class vs cell line)? Per experiment (audit_I19 inhibitor experiments binder-masked; stimulation experiments):
responders |log2FC|>1 in a direction (>=20 S/T sites) vs non-responders |log2FC|<0.3; per S/T kinase (Cantley, 303):
log OR of motif pct>=90. Profiles compared between experiment pairs by Pearson r, grouped by pair type
(same/different class x same/different cell line x same/different study)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import load_exp, mask_binders, CA as ACA
CA = '/workspace/results/b_I05/cache/'; OUT = '/workspace/results/b_I05/tables/'
MP = pd.read_parquet(CA + 'motif_pct.parquet')
S = pd.read_csv('/workspace/results/infra_compendium/sites.tsv', sep='\t', usecols=['site', 'residue'], low_memory=False).set_index('site')
KST = pd.read_csv('<HACKATHON_ROOT>/kinaseSubMain/external/kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases/kinase_data/kinome_information.tsv', sep='\t')
kst = [k for k in KST[KST.TYPE == 'ser_thr'].GENE_NAME.unique() if k in MP.columns]
prof = []; meta = []
for kind in ['inh', 'stim']:
    if kind == 'inh':
        X, E = load_exp('inh'); A = pd.DataFrame(mask_binders(X, E), index=X.index, columns=X.columns)
    else:
        M = pd.read_csv(ACA + 'stim_meta.tsv', sep='\t', dtype={'study_pmid': str})
        drop = (M.control_raw == 'LPA') | (M.study_pmid == '29020638') | M.condition_raw.str.contains('p38i', na=False) | M.agent.isin(['cytarabine', 'methotrexate'])
        M = M[~drop]; Xc = pd.read_parquet(ACA + 'stim_cond.parquet')[M.conditionID]
        M['exp'] = M.study_pmid + '|' + M.cell_line_norm.astype(str) + '|' + M.cls
        A = Xc.T.groupby(M.exp.values).mean().T
        E = M.groupby('exp').agg(study=('study_pmid', 'first'), cell=('cell_line_norm', 'first'), cls=('cls', 'first')).reindex(A.columns).reset_index()
    A = A[S.residue.reindex(A.index).isin(['S', 'T']).values]
    Mk = (MP.reindex(A.index)[kst].values >= 90)
    for j in range(A.shape[1]):
        x = A.iloc[:, j].values; fin = ~np.isnan(x)
        if fin.sum() < 500: continue
        prim_sg = -1 if kind == 'inh' else 1
        bg = fin & (np.abs(np.nan_to_num(x)) < 0.3)
        for lab, sg in [('primary', prim_sg), ('reverse', -prim_sg)]:
            rsp = fin & (sg * np.nan_to_num(x) > 1)
            if rsp.sum() < 20: continue
            a = Mk[rsp].sum(0) + .5; b = rsp.sum() - Mk[rsp].sum(0) + .5; c = Mk[bg].sum(0) + .5; d = bg.sum() - Mk[bg].sum(0) + .5
            prof.append(np.log(a * d / (b * c)))
            meta.append(dict(kind=kind, exp=E.exp.iloc[j] if 'exp' in E else A.columns[j], study=E.study.iloc[j], cell=str(E.cell.iloc[j]), cls=E.cls.iloc[j],
                             dir=lab, n_resp=int(rsp.sum())))
P = pd.DataFrame(np.array(prof), columns=kst); Mt = pd.DataFrame(meta)
P.to_parquet(CA + 'kinase_enrich_profiles.parquet'); Mt.to_csv(CA + 'kinase_enrich_meta.tsv', sep='\t', index=False)
# ---- top enriched kinases among reverse responders (study-level mean, per kind), and per class
rows = []
for (kind, dr), m in Mt.groupby(['kind', 'dir']):
    pp = P.loc[m.index]
    st = pp.groupby(m.study.values).mean()
    for k in kst:
        rows.append(dict(kind=kind, dir=dr, kinase=k, mean_lor=st[k].mean(), frac_studies_pos=(st[k] > 0).mean(), n_studies=len(st), frac_exps_pos=(pp[k] > 0).mean()))
K = pd.DataFrame(rows)
K.to_csv(OUT + 'T5_kinase_enrichment_by_direction.tsv', sep='\t', index=False)
for (kind, dr), d in K.groupby(['kind', 'dir']):
    print(kind, dr, d.n_studies.iloc[0], 'studies; top:', ', '.join(f'{r.kinase} {r.mean_lor:.2f} ({r.frac_studies_pos:.0%})' for r in d.nlargest(12, 'mean_lor').itertuples()))
# per class top reverse kinases
pcl = []
for (kind, c), m in Mt[Mt.dir == 'reverse'].groupby(['kind', 'cls']):
    if m.study.nunique() < 2: continue
    st = P.loc[m.index].groupby(m.study.values).mean().mean()
    top = st.nlargest(6)
    pcl.append(dict(kind=kind, cls=c, n_exp=len(m), n_studies=m.study.nunique(), top=', '.join(f'{k} {v:.2f}' for k, v in top.items())))
PCL = pd.DataFrame(pcl); PCL.to_csv(OUT + 'T5b_reverse_top_kinases_by_class.tsv', sep='\t', index=False)
print(PCL.to_string())
# ---- consistency of profiles by pair type
Z = P.values - P.values.mean(1, keepdims=True)  # remove global shift; correlate shapes
Z = Z / np.linalg.norm(Z, axis=1, keepdims=True)
Rm = Z @ Z.T
out = []
for kind in ['inh', 'stim']:
    for dr in ['primary', 'reverse']:
        ix = Mt.index[(Mt.kind == kind) & (Mt.dir == dr)].values
        m = Mt.loc[ix]
        cls = m.cls.values; cell = m.cell.values; st = m.study.values
        I, J = np.triu_indices(len(ix), 1)
        r = Rm[np.ix_(ix, ix)][I, J]
        typ = np.select([(cls[I] == cls[J]) & (cell[I] == cell[J]) & (st[I] != st[J]), (cls[I] == cls[J]) & (cell[I] != cell[J]) & (st[I] != st[J]),
                         (cls[I] != cls[J]) & (cell[I] == cell[J]) & (st[I] != st[J]), (cls[I] != cls[J]) & (cell[I] != cell[J]) & (st[I] != st[J]),
                         (cls[I] == cls[J]) & (cell[I] != cell[J]) & (st[I] == st[J]), (cls[I] != cls[J]) & (cell[I] == cell[J]) & (st[I] == st[J]),
                         (cls[I] != cls[J]) & (cell[I] != cell[J]) & (st[I] == st[J])],
                        ['xS sameClass sameCell', 'xS sameClass diffCell', 'xS diffClass sameCell', 'xS diffClass diffCell',
                         'wS sameClass diffCell', 'wS diffClass sameCell', 'wS diffClass diffCell'], 'other')
        df = pd.DataFrame(dict(r=r, typ=typ, si=st[I], sj=st[J]))
        for t, g in df.groupby('typ'):
            # study-pair weighted mean + study bootstrap
            spm = g.groupby(['si', 'sj']).r.mean()
            stl = np.unique(np.r_[spm.index.get_level_values(0), spm.index.get_level_values(1)])
            rng = np.random.default_rng(5); bs = []
            a_ = spm.index.get_level_values(0).values; b_ = spm.index.get_level_values(1).values
            for _ in range(500):
                w = pd.Series(rng.choice(stl, len(stl))).value_counts()
                ww = w.reindex(a_).fillna(0).values * w.reindex(b_).fillna(0).values
                if ww.sum() > 0: bs.append(np.average(spm.values, weights=ww))
            out.append(dict(kind=kind, dir=dr, pair_type=t, mean_r=spm.mean(), lo=np.percentile(bs, 2.5) if bs else np.nan, hi=np.percentile(bs, 97.5) if bs else np.nan,
                            n_pairs=len(g), n_studypairs=len(spm)))
O = pd.DataFrame(out); O.to_csv(OUT + 'T5c_profile_consistency.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(O.round(3).to_string())
