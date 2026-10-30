"""Task 4: are reverse-direction (up) dose responses to kinase inhibitors less potent than primary (down) ones?
Reuses b_I17 decryptM 4PL fits (read-only): significant curves (sig, EC50 in range, |delta|>=log2 1.5).
Per experiment and |delta| bin: median pEC50(down) - median pEC50(up) (>0 = primary more potent); also Hill.
Aggregation: experiment -> agent mean -> mean over agents, bootstrap over agents. Circularity: none (no selection on potency).
Control: per experiment, shuffle direction labels within |delta| bin (expected 0)."""
import numpy as np, pandas as pd
OUT = '/workspace/results/b_I05/tables/'
c = pd.read_parquet('/workspace/results/b_I17/curves_all.parquet', columns=['exp_id', 'site', 'delta', 'pEC50', 'hill', 'sig', 'agent', 'cell_line', 'time_min', 'pclass'])
c = c[c.sig & (c.pclass == 'kinase_inhibitor')].copy()
c['dir'] = np.where(c.delta < 0, 'down', 'up'); c['bin'] = pd.cut(c.delta.abs(), [0.58, 1, 2, 99], labels=['0.58-1', '1-2', '>2'])
rng = np.random.default_rng(6)
rows = []
for (e, b), g in c.groupby(['exp_id', 'bin'], observed=True):
    dn, up = g[g.dir == 'down'], g[g.dir == 'up']
    if len(dn) < 5 or len(up) < 5: continue
    lab = np.array(g.dir.values, dtype=object); nl = []
    for _ in range(100):
        rng.shuffle(lab); nl.append(np.median(g.pEC50.values[lab == 'down']) - np.median(g.pEC50.values[lab == 'up']))
    rows.append(dict(exp_id=e, agent=g.agent.iloc[0], cell=g.cell_line.iloc[0], time=g.time_min.iloc[0], bin=b, n_down=len(dn), n_up=len(up),
                     d_pEC50=dn.pEC50.median() - up.pEC50.median(), d_hill=dn.hill.median() - up.hill.median(), null_sd=np.std(nl)))
R = pd.DataFrame(rows); R.to_csv(OUT + 'T6_dose_potency_per_exp.tsv', sep='\t', index=False)
out = []
for b in ['0.58-1', '1-2', '>2', 'all']:
    d = R if b == 'all' else R[R.bin == b]
    for col in ['d_pEC50', 'd_hill']:
        ag = d.groupby('agent')[col].mean()
        bs = [ag.sample(len(ag), replace=True, random_state=rng.integers(1e9)).mean() for _ in range(2000)]
        out.append(dict(bin=b, metric=col, est=ag.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), n_agents=len(ag), n_exp=d.exp_id.nunique(),
                        frac_agents_pos=(ag > 0).mean()))
O = pd.DataFrame(out); O.to_csv(OUT + 'T6_dose_potency_summary.tsv', sep='\t', index=False)
print(O.round(3).to_string())
print(R[R.bin == '1-2'].sort_values('d_pEC50').round(2).to_string())
