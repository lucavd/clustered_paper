#!/bin/bash
# s8b_sub_pipeline.sh -- S8b sub-task: regenerate Table S36, Table S37, the MD/KC rows of Table S26, the MCSE caption
# values, Figure S1, the keys file tables/V5_S8b_sub_keys.csv and the snippets. Read-only on results_v4/ and every
# existing table; writes only tables/V5_S8b_sub_*.csv, reports/morris/S8b_figures/figS1.*, reports/morris/S8b_logs/sub/.
# One process at a time (shared host). Run from anywhere: bash ~/MMM/R/s8b/sub/s8b_sub_pipeline.sh
set -u
cd ~/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=reports/morris/S8b_logs/sub; mkdir -p $L reports/morris/S8b_figures
run() { b=$1; shift; t0=$(date +%s)
  if "$@" > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -20 $L/$b.log; exit 1; fi; }
run s8b_table_s36     Rscript --vanilla R/s8b/sub/s8b_table_s36.R
run s8b_table_s37     Rscript --vanilla R/s8b/sub/s8b_table_s37.R
run s8b_mdkc_s26      Rscript --vanilla R/s8b/sub/s8b_mdkc_s26.R
run s8b_mcse_captions python3 -B py/s8b/s8b_mcse_captions.py
run s8b_figS1         python3 -B py/s8b/s8b_figS1.py
run s8b_sub_keys      Rscript --vanilla R/s8b/sub/s8b_sub_keys.R
if [ -f py/s8b/s8b_sub_snippets.py ]; then run s8b_sub_snippets python3 -B py/s8b/s8b_sub_snippets.py; fi
md5sum tables/V5_S8b_sub_*.csv reports/morris/S8b_figures/figS1.png > $L/md5_outputs.txt
echo "sub pipeline done"
