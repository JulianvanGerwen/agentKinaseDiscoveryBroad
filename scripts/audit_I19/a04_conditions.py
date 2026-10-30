"""(c) time/dose matching and (e) noise ceiling / context decomposition at the condition level (no averaging over time/dose).
Inhibitor conditions, raw vs control, own-binder proteins masked; sites measured in >=6 conditions."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/scripts/audit_I19'); from common_audit import *
OUT = '/workspace/results/audit_I19/'
pd.set_option('display.width', 250)
X, E = load_cond('inh'); keep = X.notna().sum(axis=1).values >= 6; X = X[keep]
A = mask_binders(X, E)
st, cl, ag, ce = E.study.values, E.cls.values, E.agent.values, E.cell.astype(str).values
t = E.time_min.values; dose = E.dose_nM.values
tb = np.select([t <= 30, (t > 30) & (t < 120), t >= 120], ['early', 'mid', 'late'], 'NA'); tb[np.isnan(t)] = 'NA'
allcls = {c: (E.cls == c).values for c in E.cls.unique()}
multi = {c: v for c, v in allcls.items() if E[E.cls == c].study.nunique() >= 2}
rows = []
def run(label, R):
    if len(R) == 0: print(label, 'no pairs'); return
    S = summarise(R)
    out = S[['dir', 'bin', 'pairs', 'study_pairs', 'rep', 'base', 'opp', 'opp_base', 'lor', 'meanB']].copy(); out.insert(0, 'test', label)
    # within-study tests: study-pair weighting collapses to per-study; also report directional excess
    out['dir_excess'] = (out.rep - out.base) - (out.opp - out.opp_base)
    for b in ['1-2', '2-99']:
        try:
            o, lo, hi, ns, nsp = boot_diff(R, 'lor', b); out.loc[out.bin == b, 'dLOR'] = o; out.loc[out.bin == b, 'dLOR_CI'] = f'[{lo:.2f},{hi:.2f}]'
        except Exception: pass
    rows.append(out); print(out.round(3).to_string(), flush=True)
# (c) time-matched cross-study pairs
for b in ['early', 'mid', 'late']:
    run(f'c_time_both_{b}', pair_rows(A, E, multi, null_pool=None, pair_filter=lambda i, j, b=b: tb[i] == b and tb[j] == b))
run('c_time_mismatched', pair_rows(A, E, multi, null_pool=None, pair_filter=lambda i, j: tb[i] != tb[j] and 'NA' not in (tb[i], tb[j])))
run('c_dose_both_>=1uM', pair_rows(A, E, multi, null_pool=None, pair_filter=lambda i, j: dose[i] >= 1000 and dose[j] >= 1000))
run('c_dose_both_<1uM_or_one', pair_rows(A, E, multi, null_pool=None, pair_filter=lambda i, j: dose[i] < 1000 or dose[j] < 1000))
run('c_time_both_<=60_dose_both_>=1uM', pair_rows(A, E, multi, null_pool=None, pair_filter=lambda i, j: t[i] <= 60 and t[j] <= 60 and dose[i] >= 1000 and dose[j] >= 1000))
# (e) noise ceiling & decomposition, within-study pairs (all classes), plus cross-study same class
W = dict(cross_study=False, null_pool=None)
run('e1_same_study_same_drug_same_cell_diff_dose_or_time', pair_rows(A, E, allcls, pair_filter=lambda i, j: st[i] == st[j] and ag[i] == ag[j] and ce[i] == ce[j], **W))
run('e2_same_study_same_class_diff_drug_same_cell', pair_rows(A, E, allcls, pair_filter=lambda i, j: st[i] == st[j] and ag[i] != ag[j] and ce[i] == ce[j], **W))
run('e3_same_study_same_drug_diff_cell', pair_rows(A, E, allcls, pair_filter=lambda i, j: st[i] == st[j] and ag[i] == ag[j] and ce[i] != ce[j], **W))
run('e4_same_study_same_class_diff_drug_diff_cell', pair_rows(A, E, allcls, pair_filter=lambda i, j: st[i] == st[j] and ag[i] != ag[j] and ce[i] != ce[j], **W))
run('e5_cross_study_same_class_same_cell', pair_rows(A, E, allcls, pair_filter=lambda i, j: st[i] != st[j] and ce[i] == ce[j], **W))
run('e6_cross_study_same_class_diff_cell', pair_rows(A, E, allcls, pair_filter=lambda i, j: st[i] != st[j] and ce[i] != ce[j], **W))
# null: same study, same cell, different class (pairs across classes -> one pseudo-group of all exps)
allg = {'ALL': np.ones(len(E), bool)}
rng = np.random.default_rng(1)
sub = np.zeros(len(E), bool)
run('e0_NULL_same_study_same_cell_diff_class', pair_rows(A, E, allg, pair_filter=lambda i, j: st[i] == st[j] and ce[i] == ce[j] and cl[i] != cl[j] and (i * 7 + j) % 5 == 0, **W))
run('e0b_NULL_cross_study_diff_class_same_cell', pair_rows(A, E, allg, pair_filter=lambda i, j: st[i] != st[j] and ce[i] == ce[j] and cl[i] != cl[j], **W))
pd.concat(rows).to_csv(OUT + 'T_tests_c_e.tsv', sep='\t', index=False)
