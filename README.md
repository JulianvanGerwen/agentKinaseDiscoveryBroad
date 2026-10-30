# Agentic AI for bioinformatics: searching for new kinase-signalling biology in phosphoproteomics data

> **This repository is an exercise in testing agentic AI for bioinformatics.** A team of AI agents ran the project end to end from a short human-written brief: a lead "PI" agent (Claude Opus 5.5) directing Claude Opus 5.5 and Sonnet 5.5 sub-agents. The agents wrote all analyses, code, figures and text. **No result has been verified by a human.** Treat every biological claim as a hypothesis, and see the limitations below.

<!-- Placeholder for a human-written introduction (context of the exercise, verification status). The agents did not write anything here. -->

## Start here: the research briefing

### **[Read the briefing (PDF)](docs/briefing/BRIEFING.pdf)** · [Markdown version](docs/briefing/BRIEFING.md)

The briefing is the main way to engage with this project. It is a paper-length account for a broad scientific audience. It covers the project both as **an exercise in agentic science** (how the AI team was organised, how claims were checked, what went wrong) and as **a bioinformatics study** (what was found, how confident we are, and what was tried and dropped). It explains its terminology, credits the literature ideas it builds on, and includes eight figures.

## The project in one paragraph

The agents analysed:
- a compendium of 1,937 phosphoproteomics perturbation experiments from 178 studies;
- CPTAC tumour and normal-tissue phosphoproteomes from 10 cancer types;
- several public datasets.

They used two kinase-analysis tools built by earlier AI-run projects. The goal was "several novel and compelling discoveries about kinase signalling". Over about 18 hours the PI agent ran 18 discovery branches, 4 targeted validations, 6 independent audits and 3 reviews. The outcome, by audited confidence:

| Tier | Finding |
|---|---|
| **1: discovery** | Phosphorylation of **KPNA4 S60** (importin-α3) depends on **CK2**, not on the MAPK pathway as previously published. Our compendium data are supported by an independent public chemical-genetic rescue. Direct phosphorylation is not yet shown. |
| 2: qualified finding | Signalling programmes learned from cell-culture perturbations reappear as coherent units in human tumours *and* normal tissue; tumours mainly gain a more coordinated DNA-damage programme. |
| 2: qualified finding | Phosphorylation changes in a perturbation's main direction reproduce across laboratories; reverse-direction changes mostly do not (a methodological caution). |
| 2: qualified finding | In tumours, receptor signalling driven by receptor overexpression is buffered before ERK; cancer lesions shift pathway output levels but not pathway sensitivity. |
| Caveated | Candidate mTORC1 sites on its lysosomal docking machinery (LAMTOR1, RRAGC); a phosphosite "turnover" property; a drug-potency/ATP consistency check; an audit of curated kinase–substrate annotations. |
| Negative / dropped | Nine of 18 discovery branches, e.g. understudied kinases, phosphatase lesions, copy-number causal tests, allele-specific signalling, an ERK "ceiling" in EGFR-mutant lung cancer that failed replication. |

Overall, the brief was **only partly met**: one solid but modest named discovery and three qualified findings, from a broad and largely negative search.

## Guide to the repository

**Main documents**

| File | What it is |
|---|---|
| [`docs/briefing/`](docs/briefing/) | **The briefing** (PDF, Markdown, HTML), its figures (`figures/`, with plotting scripts and `data_sources.md`), and its independent fact-check (`FACTCHECK.md`) |
| [`PROMPT.md`](PROMPT.md) | The human-written brief, verbatim |
| [`REPORT.md`](REPORT.md) | The technical report: ranked findings, evidence, limitations |
| [`RESEARCH_TREE.md`](RESEARCH_TREE.md) | The PI agent's research record: every branch, decision (D0–D9) and milestone summary |
| [`LESSONS.md`](LESSONS.md) | Process rules the PI agent kept, each with the evidence that motivated it |
| [`LITERATURE.md`](LITERATURE.md), [`docs/lit/`](docs/lit/) | Literature review and the 97-row novelty checklist used to judge novelty |
| [`docs/IDEAS.md`](docs/IDEAS.md) | The idea pool and its ranking |
| [`docs/retro/`](docs/retro/) | Milestone retrospectives by fresh agents, and the final fact-check of `REPORT.md` |
| [`DATA.md`](DATA.md) | Data register for the compendium and CPTAC assets built by the agents |

**Analysis outputs and code**

| Path | What it holds |
|---|---|
| `results/b_*/SUMMARY.md` | One summary per discovery branch (e.g. `b_I03` new kinase–substrate edges, `b_Q8` tumour vs normal tissue) |
| `results/audit_*/AUDIT.md` | Independent audits that recomputed and challenged each candidate finding |
| `results/val_*/SUMMARY.md` | Targeted validations (cross-validation, measured ischaemia times, KPNA4/mTORC1 tests) |
| `results/infra_*` | Summaries of the compendium and CPTAC infrastructure |
| `scripts/<branch>/` | Analysis code, one folder per branch |

## What is not included

- **Raw data are not included.** This covers the phosphoproteomics compendium, PhosphoSitePlus, CPTAC/PDC/GDC, DepMap, GDSC, Kinobeads, PRIDE PXD038050 and the other datasets; their sources and licences are in [`NOTICE.md`](NOTICE.md).
- **Patient-level tables and large matrices** (`*.parquet`, `cache/`) are excluded by `.gitignore`.
- **The two tool repositories the agents used** (the kinase–substrate predictor *kinaseSubMain* and the activity estimator *kinaseActivityMain*) are separate projects. Their location is shown as the placeholder `<HACKATHON_ROOT>`. See [`docs/PATHS.md`](docs/PATHS.md) for the placeholders and what a reader needs to set to re-run scripts.

To rebuild the briefing PDF:

```bash
pip install markdown xhtml2pdf
python docs/briefing/render_pdf.py
```

## Limitations to read before using anything here

1. **No human verification yet.** The agents audited each other, and audits changed the wording of every positive claim, but independent agent audits are not a substitute for expert review.
2. **The tools have known limits.** Kinase activity inference works at the level of kinase *families*, and tyrosine phosphorylation is sparse in tumour data.
3. **The inhibitor data are concentrated.** They come largely from two laboratories and from short treatments of cultured cells.
4. **Tumour analyses are cross-sectional.** Pre-analytical effects beyond recorded cold-ischaemia time cannot be fully excluded.
5. **The KPNA4 finding needs biochemical confirmation.** It needs an in vitro CK2 kinase assay and a short time course with the inhibitor-resistant CK2 rescue.

## Licences

| Material | Licence |
|---|---|
| Code (`scripts/`, `docs/briefing/*.py`, figure scripts) | MIT; see [`LICENSE`](LICENSE) |
| Documentation and briefing | CC BY 4.0 |
| Derived tables with PhosphoSitePlus-derived kinase–substrate information | CC BY-NC-SA 3.0, with attribution to PhosphoSitePlus (www.phosphosite.org) |

Third-party data are not redistributed; see [`NOTICE.md`](NOTICE.md).
