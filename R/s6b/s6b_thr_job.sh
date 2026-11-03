#!/bin/bash
# S6b part A run: parity test (GA1), md5 before (GA2), fits for threshold 0.01 (three DGPs, 48 cores, one thread each).
set -u
cd ~/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results_v4/S6b reports/morris/S6b_logs
L=reports/morris/S6b_logs
Rscript --vanilla R/testing/test_h0h1_parity.R > $L/GA1_parity.log 2>&1; echo "GA1 parity exit $?"; tail -3 $L/GA1_parity.log
( cd results_full && find . -name '*.rds' -type f | sort | xargs md5sum ) > results_v4/S6b/md5_results_full_before.txt
md5sum results_v4/S6cc/cc_all.rds results_v4/variant_cc/S4/rows_v4.rds > results_v4/S6b/md5_inputs_before.txt
wc -l < results_v4/S6b/md5_results_full_before.txt
for D in GenBinaryY GLMM-native GLMM-AR1; do Rscript --vanilla R/s6b/s6b_thr_run.R $D 48 2>&1 | grep -v 'Cgee\|running glm' | tail -3; done
