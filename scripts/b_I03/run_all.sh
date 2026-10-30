#!/bin/bash
# I03 pipeline (~25 min, <=4 threads, <20 GB). Env: /home/plm_env/bin/python
set -e; P=/home/plm_env/bin/python; D=/workspace/scripts/b_I03
for s in 01_pert 02_dose 03_tumour 04_context 04b_kinasesub_full 05_assemble 06_validate 06b_interaction_model 06c_conditional_or 07_nominate 08_edge_detail 09_figures_tables; do
  $P $D/$s.py; done
