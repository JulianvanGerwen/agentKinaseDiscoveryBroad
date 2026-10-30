"""Task 2/3: reduced-form (CNV -> protein-adjusted site), Wald/2SLS ratio, naive correlation, controls.
Per cohort, instrumented 'genes' = strong kinases (F>10, b>0) + up to 400 strong non-kinase fake genes.
Outputs:
  cache/<C>/rf.npz          t_rf, beta_rf, t_naive, t_rf_focal, t_rf_nbadj for kinases (K x S), sites, kinases
  tables/T2_eval.tsv        per (gene, cohort, instrument-type): AUROC vs PSP / motif (own), specificity ranks
  tables/T2_null.tsv        null AUROCs (fake genes, permutations)
"""
import numpy as np, pandas as pd, sys
from scipy import stats
from lib import *
rng = np.random.default_rng(0)
NPERM = 20

fs = pd.read_parquet(f'{RES}/tables/T1_first_stage_all_genes.parquet')
ann = pd.read_parquet(f'{CACHE}/gene_annot.parquet').drop_duplicates('uniprot').set_index('uniprot')
psp = pd.read_parquet(f'{CACHE}/psp_human.parquet')
psp['kin'] = psp.KIN_ACC_ID.str.split('-').str[0]
psp['site_key'] = psp.SUB_ACC_ID.str.split('-').str[0] + '_' + psp.SUB_MOD_RSD
motif = pd.read_parquet(f'{CACHE}/motif_pct.parquet')
sym = ann.symbol.to_dict()
u_of_sym = {v: k for k, v in sym.items()}
# PSP kinase gene set dictionary (by uniprot)
psp_sets = psp.groupby('kin').site_key.apply(set).to_dict()
psp_kin_list = [k for k, v in psp_sets.items() if len(v) >= 5]
motif_kin_list = list(motif.columns)

evals, nulls = [], []
for C in (sys.argv[1:] or COHORTS):
    padj = load(C, 'padj'); prot = load(C, 'prot'); cnv = load(C, 'cnv'); foc = load(C, 'cnv_focal'); cov = load(C, 'cov')
    f = fs[fs.cohort == C].set_index('uniprot')
    strongK = f[(f.is_kinase) & (f.F_cnv > 10) & (f.b_cnv > 0)].index.intersection(cnv.index)
    fakeall = f[(~f.is_kinase) & (f.F_cnv > 10) & (f.b_cnv > 0)].index.intersection(cnv.index)
    fake = rng.choice(fakeall, min(400, len(fakeall)), replace=False)
    genes = list(strongK) + list(fake)
    cv = cov.values
    Z = resid_cov(cnv.loc[genes].T.values.astype(float), cv)            # n x G
    Zf = resid_cov(foc.loc[strongK].T.values.astype(float), cv)
    Pk = prot.loc[strongK].T.values.astype(float)
    Pk = np.where(np.isfinite(Pk), Pk, np.nanmean(Pk, 0))
    Pk = resid_cov(Pk, cv)
    Y = padj.T.values.astype(float)                                      # n x S
    sites = padj.index.values
    site_u = np.array([s.split('_')[0] for s in sites])
    site_arm = np.array(pd.Series(site_u).map(ann.arm).astype(object).values)
    b_rf, t_rf, n_s = masked_t(Z, Y)
    _, t_nv, _ = masked_t(Pk, Y)
    _, t_fo, _ = masked_t(Zf, Y)
    # first-stage slope on residualised Z (for Wald ratio)
    Pall = prot.reindex(genes).T.values.astype(float)
    bfs = np.array([np.polyfit(Z[np.isfinite(Pall[:, j]), j], Pall[np.isfinite(Pall[:, j]), j], 1)[0] for j in range(len(genes))])
    # neighbour-adjusted reduced form: partial out the mean standardised protein of co-amplified non-kinase neighbours
    t_nb = np.full((len(strongK), len(sites)), np.nan); nb_n = []
    cnv_c = cnv.loc[cnv.index.intersection(prot.index)]
    cnv_arm = cnv_c.index.map(ann.arm)
    for i, k in enumerate(strongK):
        a = ann.arm.get(k)
        cand = cnv_c[(cnv_arm == a) & (cnv_c.index != k) & ~cnv_c.index.isin(f.index[f.is_kinase])]
        if len(cand) == 0:
            nb_n.append(0); continue
        r = np.corrcoef(cnv.loc[k].values, cand.values)[0, 1:]
        nb = cand.index[r > 0.8]
        nb = [g for g in nb if g in f.index and f.loc[g, 'F_cnv'] > 10]   # neighbours whose protein follows the CNV
        nb_n.append(len(nb))
        if len(nb) == 0:
            continue
        W = prot.loc[nb].T
        W = ((W - W.mean()) / W.std()).mean(1).fillna(0).values
        W = resid_cov(W[:, None], cv)
        X = np.column_stack([Z[:, i], W[:, 0]])
        # residualise both Z_k and Y on W (masked FWL approximation: Y via simple masked regression)
        zr = Z[:, i] - W[:, 0] * (W[:, 0] @ Z[:, i]) / (W[:, 0] @ W[:, 0])
        bw, _, _ = masked_t(W, Y)
        Yr = Y - W @ bw
        _, tt, _ = masked_t(zr[:, None], Yr)
        t_nb[i] = tt[0]
    K = len(strongK)
    np.savez_compressed(f'{CACHE}/{C}/rf.npz', t_rf=t_rf[:K].astype(np.float32), b_rf=b_rf[:K].astype(np.float32),
                        t_nv=t_nv.astype(np.float32), t_fo=t_fo.astype(np.float32), t_nb=t_nb.astype(np.float32),
                        n_s=n_s, sites=sites, kin=np.array(strongK), bfs=bfs[:K], nb_n=np.array(nb_n))
    # ----- evaluation -----
    # label matrices over sites: PSP (kinases with >=5 subs) and motif>=90
    Lp = np.array([[s in psp_sets[k] for s in sites] for k in psp_kin_list]).T.astype(float)  # S x Kp
    mot = motif.reindex(sites)
    Lm = (mot.values >= 90).astype(float); Vm = np.isfinite(mot.values)       # residue-compatible
    def auc_matrix(T, L, valid=None, excl=None):
        """T: (G x S) scores; L: (S x J) labels; returns G x J AUROC; excl: (G x S) bool sites to drop"""
        T = T.copy()
        if excl is not None:
            T[excl] = np.nan
        out = np.full((T.shape[0], L.shape[1]), np.nan)
        for g in range(T.shape[0]):
            ok = np.isfinite(T[g])
            r = stats.rankdata(T[g, ok])
            Lg = L[ok]
            if valid is not None:
                # per-label validity: use only residue-compatible sites -> rank within valid subset
                Vg = valid[ok]
                for j in range(L.shape[1]):
                    v = Vg[:, j]; rr = stats.rankdata(T[g, ok][v]); l = Lg[v, j] > 0
                    npos, nneg = l.sum(), (~l).sum()
                    if npos >= 3 and nneg > 0:
                        out[g, j] = (rr[l].sum() - npos * (npos + 1) / 2) / (npos * nneg)
                continue
            npos = Lg.sum(0); nneg = ok.sum() - npos
            with np.errstate(invalid='ignore', divide='ignore'):
                out[g] = np.where(npos >= 3, (r @ Lg - npos * (npos + 1) / 2) / (npos * nneg), np.nan)
        return out
    # exclusion masks: own-protein sites and same-arm sites
    garm = np.array([ann.arm.get(g) for g in genes], dtype=object)
    excl = (site_u[None, :] == np.array(genes)[:, None]) | (site_arm[None, :] == garm[:, None])
    exclK = excl[:K]
    A_psp = auc_matrix(t_rf, Lp, excl=excl)                 # (K+fake) x Kp
    A_psp_nv = auc_matrix(t_nv, Lp, excl=exclK)
    A_psp_fo = auc_matrix(t_fo, Lp, excl=exclK)
    A_psp_nb = auc_matrix(t_nb, Lp, excl=exclK)
    # motif: residue class handled by NaN in motif -> treat as negatives only within compatible: approximate by
    # restricting to compatible residue via per-kinase split (ST vs Y)
    isY = np.array([s.split('_')[1][0] == 'Y' for s in sites])
    kY = np.isfinite(mot.values[isY]).any(0) if isY.any() else np.zeros(mot.shape[1], bool)
    kY = ~np.isfinite(mot.values[~isY]).any(0)   # kinase scored only on Y sites -> Tyr kinase
    def motif_auc(T, ex):
        res = np.full((T.shape[0], Lm.shape[1]), np.nan)
        for cls, kmask in [(isY, kY), (~isY, ~kY)]:
            ex2 = ex | ~cls[None, :]
            a = auc_matrix(T, np.nan_to_num(Lm) * cls[:, None], excl=ex2)
            res[:, kmask] = a[:, kmask]
        return res
    A_mot = motif_auc(t_rf, excl); A_mot_nv = motif_auc(t_nv, exclK); A_mot_fo = motif_auc(t_fo, exclK); A_mot_nb = motif_auc(t_nb, exclK)
    # permutation null for kinases: permute samples of Z
    A_perm_psp = []; A_perm_mot = []
    for p in range(NPERM):
        Zp = Z[rng.permutation(len(Z)), :K]
        _, tp, _ = masked_t(Zp, Y)
        A_perm_psp.append(auc_matrix(tp, Lp, excl=exclK)); A_perm_mot.append(motif_auc(tp, exclK))
    pk_idx = {k: j for j, k in enumerate(psp_kin_list)}
    mk_idx = {k: j for j, k in enumerate(motif_kin_list)}
    gene_chr = np.array([ann.chrom.get(g) for g in genes])
    for i, g in enumerate(genes):
        is_k = i < K
        s = sym.get(g, g)
        rec = dict(cohort=C, uniprot=g, gene=s, is_kinase=is_k, F=f.loc[g, 'F_cnv'], F_focal=f.loc[g, 'F_focal'],
                   b_fs=bfs[i], n_sig_rf=int((bh(2 * stats.t.sf(np.abs(t_rf[i]), n_s - 3)) < 0.1).sum()))
        if is_k:
            rec['nb_n'] = nb_n[i]
            jp = pk_idx.get(g); jm = mk_idx.get(s)
            if jp is not None:
                rec.update(auc_psp_rf=A_psp[i, jp], auc_psp_nv=A_psp_nv[i, jp], auc_psp_fo=A_psp_fo[i, jp], auc_psp_nb=A_psp_nb[i, jp],
                           n_psp=int((Lp[:, jp] * ~exclK[i]).sum()),
                           psp_spec_rank=stats.percentileofscore(A_psp[i][np.isfinite(A_psp[i])], A_psp[i, jp]),
                           psp_spec_rank_nv=stats.percentileofscore(A_psp_nv[i][np.isfinite(A_psp_nv[i])], A_psp_nv[i, jp]),
                           perm_psp=np.nanmean([a[i, jp] for a in A_perm_psp]), perm_psp_sd=np.nanstd([a[i, jp] for a in A_perm_psp]))
                fk = np.where((np.arange(len(genes)) >= K) & (gene_chr != gene_chr[i]))[0]
                rec['fake_psp_mean'] = np.nanmean(A_psp[fk, jp]); rec['fake_psp_p'] = (np.sum(A_psp[fk, jp] >= A_psp[i, jp]) + 1) / (np.isfinite(A_psp[fk, jp]).sum() + 1)
                rec['fake_psp_mean_nv'] = np.nan
            if jm is not None:
                rec.update(auc_mot_rf=A_mot[i, jm], auc_mot_nv=A_mot_nv[i, jm], auc_mot_fo=A_mot_fo[i, jm], auc_mot_nb=A_mot_nb[i, jm],
                           mot_spec_rank=stats.percentileofscore(A_mot[i][np.isfinite(A_mot[i])], A_mot[i, jm]),
                           mot_spec_rank_nv=stats.percentileofscore(A_mot_nv[i][np.isfinite(A_mot_nv[i])], A_mot_nv[i, jm]),
                           perm_mot=np.nanmean([a[i, jm] for a in A_perm_mot]), perm_mot_sd=np.nanstd([a[i, jm] for a in A_perm_mot]))
                fk = np.where((np.arange(len(genes)) >= K) & (gene_chr != gene_chr[i]))[0]
                rec['fake_mot_mean'] = np.nanmean(A_mot[fk, jm]); rec['fake_mot_p'] = (np.sum(A_mot[fk, jm] >= A_mot[i, jm]) + 1) / (np.isfinite(A_mot[fk, jm]).sum() + 1)
        evals.append(rec)
    # store full null distributions of fake genes (all kinase labels) as summary
    nulls.append(pd.DataFrame({'cohort': C, 'fake_psp_all': np.nanmean(A_psp[K:], 1), 'fake_mot_all': np.nanmean(A_mot[K:], 1)}))
    print(C, 'K', K, 'fake', len(fake), 'sites', len(sites), flush=True)

ev = pd.DataFrame(evals)
tag = '' if len(sys.argv) == 1 else '_' + '_'.join(sys.argv[1:])
ev.to_csv(f'{RES}/tables/T2_eval{tag}.tsv', sep='\t', index=False)
pd.concat(nulls).to_csv(f'{RES}/tables/T2_null_fake{tag}.tsv', sep='\t', index=False)
