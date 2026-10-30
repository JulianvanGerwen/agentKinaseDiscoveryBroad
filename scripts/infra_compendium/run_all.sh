#!/usr/bin/env bash
# Rebuild all infra_compendium assets (~2 min total, <8 GB RAM, 8 threads). Env: /home/plm_env/bin/python
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/home/plm_env/bin/python
mkdir -p logs
$PY 01_build_matrix.py      > logs/01_build_matrix.log 2>&1
$PY 02_build_edges.py       > logs/02_build_edges.log 2>&1
$PY 03_score_corev3.py      > logs/03_score_corev3.log 2>&1
$PY 04_annotate_conditions.py > logs/04_annotate.log 2>&1
$PY 05_activity_matrices.py > logs/05_activity.log 2>&1
$PY 06_sanity_check.py      > logs/06_sanity.log 2>&1
$PY 07_kinase_annotation.py > logs/07_kinase.log 2>&1
