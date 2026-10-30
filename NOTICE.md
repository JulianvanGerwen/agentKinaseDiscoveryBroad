# Licence and third-party notice

This repository is an AI-agent-run computational biology project on kinase signalling, released as an exercise in testing agentic AI for bioinformatics.

## Licences by content type

| Content | Licence |
|---|---|
| Code (`scripts/`, `*.py`, `*.sh`, `*.R`) | MIT, see `LICENSE` |
| Documentation and reports (`*.md`, `docs/`, figures) | CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) |
| Derived tables that contain PhosphoSitePlus-derived kinase-substrate information (for example kinase-substrate edge tables, site annotations and any result table keyed to PSP kinase-substrate pairs) | CC BY-NC-SA 3.0 (https://creativecommons.org/licenses/by-nc-sa/3.0/), with attribution to PhosphoSitePlus (www.phosphosite.org) |

Other small result tables are released under CC BY 4.0, except where they carry PhosphoSitePlus-derived content as described above. If in doubt, treat a table that mentions PSP kinase-substrate relationships as CC BY-NC-SA 3.0 (non-commercial use only).

Attribution for PhosphoSitePlus: Hornbeck PV et al., PhosphoSitePlus, 2014 and later releases, www.phosphosite.org. PhosphoSitePlus is a registered trademark of Cell Signaling Technology, Inc.

## Third-party data

No raw data are redistributed in this repository. The provided compendium, PhosphoSitePlus files, downloaded files and all intermediate matrices are excluded (see `.gitignore`). Patient-level tables (CPTAC case and sample identifiers) are also excluded. Readers must obtain the data from the original sources and under their terms:

- CPTAC / Proteomic Data Commons (PDC)
- NCI Genomic Data Commons (GDC)
- DepMap
- GDSC (Genomics of Drug Sensitivity in Cancer)
- Klaeger et al. 2017 (kinobeads drug-target profiling)
- Hijazi et al. 2020 (inhibitor-target data)
- decryptM
- Menyhart et al. 2023 (PRIDE PXD038050)
- Chen et al. 2020 (CPTAC LUAD, Cell)
- PhosphoSitePlus (Cell Signaling Technology), subject to its own licence

Aggregated, derived statistics computed from these sources appear in the result tables. Check the original terms before reuse.

## Fonts
The briefing PDF embeds the DejaVu Sans fonts (`docs/briefing/fonts/`), which are free software under the Bitstream Vera / DejaVu licence (https://dejavu-fonts.github.io/License.html).
