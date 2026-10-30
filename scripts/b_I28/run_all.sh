#!/bin/bash
# I28 phosphosite turnover atlas. ~4 min, <6 GB, 4 threads. Needs cache/b_I26 (scripts/b_I26/00_prep.py).
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/home/plm_env/bin/python; cd "$(dirname "$0")"; set -e
for s in 01_build 02_atlas 03_features 04_predictors 05_activity 06_cptac 07_figures; do $PY $s.py; done
