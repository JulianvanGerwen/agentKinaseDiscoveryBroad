# LESSONS — process rules (one line each, with evidence). Pruned at M2 (T+~14h) per docs/retro/M2_retro.md.

1. Every claim needs a matched null AND an independent auditor who recomputes it with own code, assigns the C score, and cross-checks against other branches' positive controls; the producing side's own validation is not independent (audits reworded A, B, D, E; audit_D resolved the KPNA4/CK2 conflict that V-I03 missed; I20 IV ≈ fake-instrument null).
2. Demand a discriminating "mirror"/orthogonal test in every audit (stimulation mirror test turned I19 into a general rule; ischaemia footprint exposed the MAPK-axis confound in B).
3. Cross-validate single-dataset candidates in independent data as the FIRST step of a follow-up, before building on them (V1: 3/4 failed in ~7 min).
4. Count independent LABS, not studies, as replication; check exact site localisation across datasets (LAMTOR1 "2 studies" = 1 lab; S97 vs S98 conflated).
5. Test pre-analytical confounds with MEASURED variables (GDC excision-to-freezing time), not with signature proxies that may be co-readouts of the biology (the compendium 'ischaemia footprint' tracked the MAPK axis, not measured time: ρ −.03, n=723 — using it over-adjusted B by ~50%). Use activity/activating sites, not abundance, as tumour exposures (I08/I20).
6. A null from an inhibitor is only informative if on-target inhibition is verified in the same experiment (I26/audit_D: CX-4945 at 1 h doesn't move known CK2 sites in Hijazi).
7. Any "not X-driven" statement needs a script that tests X (B's "not ischaemia-driven" had no test and was wrong).
8. Hard cap: after a finding is audited, at most ONE consolidation branch; keep ≥2 unrelated branches running (Finding A absorbed I22, I23, I25 after M1 — all negative/reframing).
9. Powered nulls get one tree line and no follow-up branch (I13, I20, I23, I25, I27 negatives).
10. Sub-agents may be blocked from writing SUMMARY/AUDIT .md: require full content in the final message and save it as PI.
11. Ops: check data-access restrictions within 5 min and escalate (decryptM 2.0 403); set OMP/MKL/OPENBLAS_NUM_THREADS before python in every prompt (I23 ran ~27 cores briefly); start branches the moment their infrastructure exists (≤4 agents × 4 threads).
12. Mask a kinase's own/paralog sites before relating its sites to its activity (stored scores include autophosphorylation edges; V1).
13. Freeze the claim list and audited wording at a fixed point before write-up; the report quotes audits, not branch self-scores.
14. Report only auditor scores, and label any score that isn't one; keep a single ranking threshold across tree and report (final fact-check found PI/self scores presented as audited and a C≥4 vs C≥3.5 mismatch).
15. Write the final report from AUDIT.md wording, not tree summaries, and fact-check it independently before calling it final (fact-check found ~15 overstatements, 0 fabricated numbers).
