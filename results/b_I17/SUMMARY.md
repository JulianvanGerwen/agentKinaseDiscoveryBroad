# I17/I18/I21 SUMMARY — dose-resolved pharmacology (revised; saved by PI from agent message)
Prior work: decryptM (Zecha 2023) overlaid pKd/pEC50 only; decryptM 2.0 (Bayer 2025 bioRxiv) does cross-drug kinase→site "potency coherence" (I18 scooped); ATP shift of cellular potency known from target engagement (NanoBRET Vasta 2018; KiNativ Patricelli 2011) — not from downstream phosphosites proteome-wide.
Pipeline: 74 dose series, 4PL, permutation FDR; neg controls 0–2 hits; replicate ρ .75. 21 exps/12 drugs with affinity; 10,754 significant curves.
1. **Downstream phosphosite potency shift fits ATP competition (Cheng–Prusoff).** Target-substrate median shift vs Davis 2011 (ATP-free) Kd: all 5 ATP-site drugs 1.9–2.8 log below Kd, all inside the 1–5 mM ATP band using Patricelli native ATP affinities (EGFR 18, SRC 9.8, ABL 5.2 µM); **implied cellular [ATP] 1.5–3.6 mM (median 2.3; LODO 2.0–2.3)**; Km-aware model cuts spread 2.6× vs constant offset (permutation p=.05). Kinobeads lysate Kd is a poor ATP-free reference (~1 log off). Allosteric MEKi overshoot (+0.4 to +0.8; non-ATP minus ATP-site +1.55 log, p=.024); replicated with trametinib in 3 KRAS-mutant lines (independent Kuster dataset, Zenodo 16096119). MK2206 exception. C3–4/N3/P4 — underpowered (5 drugs, 3 Km values).
2. Motif-predicted (not curated) target substrates more potent (+0.33 log, z 7.9, LODO .27–.36). C3/N3/P3.
3. K562 BCR-ABL: downstream sites respond 16×/4× earlier than ABL autosites (not general, 3/6 exps). C3/N2/P3.
4. Dasatinib engages p38α/MK2 in cells (~1 µM; OR 2.5–4.1, 4 series). C3/N2/P3.
5. Hill ≥2 in 40% (vs 23% null), poorly reproducible; no depth dependence; down-sites steeper than up (p=.005). C2/N3/P2.
Negatives: type I vs II untestable; Km dependence across drugs underpowered (ρ .45 ns); cross-drug regression ≈ shuffle; reach restates known; decryptM 2.0 restricted (Zenodo 17533475 403).
Files: results/b_I17/ (+ atp_*.tsv, kras_replication*.tsv, fig6_atp.png); scripts/b_I17/00–14.
