#!/bin/bash
# I26 PSP perturbation audit. ~3 min, <8 GB, 4 threads.
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/home/plm_env/bin/python; cd "$(dirname "$0")"
set -e
for s in 00_prep 01_audit 02_controls 03_relations 04_alternatives 05_activity 06_activity_summary 07_figures 08_curated_corrections; do $PY $s.py; done
