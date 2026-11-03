#!/bin/bash
# s6b_thr_pipeline.sh -- S6b part A: build of the variants t4 and t2, then the S4-S6 chain of R/s6/s6cc_pipeline.sh
# on each (same scripts, MMM_VARIANT=<tag>; the M1c-only and diagnostic scripts are not needed for the threshold check).
set -u
cd ~/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
L0=reports/morris/S6b_logs; mkdir -p $L0
Rscript --vanilla R/s6b/s6b_thr_build.R > $L0/s6b_thr_build.log 2>&1 || { echo "FAIL build"; tail -20 $L0/s6b_thr_build.log; exit 1; }
tail -16 $L0/s6b_thr_build.log
for TAG in t4 t2; do
  export MMM_VARIANT=$TAG
  L=reports/morris/variant_$TAG/logs; mkdir -p $L
  for s in R/postrun/agg_v4.R R/s5/s5_build.R R/s5/s5_run.R R/s5/s5_checks.R R/s6/s6_run.R R/s6/s6_run_boby.R R/s6/s6_checks.R \
           R/s6/s6_tables.R R/s6/s6_figdata.R R/s6/s6_dp_rates.R R/s6/s6_rule.R R/s6/s6_rule_safety.R; do
    b=$(basename $s .R); t0=$(date +%s)
    if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $TAG $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $TAG $b"; tail -8 $L/$b.log; exit 1; fi
  done
done
unset MMM_VARIANT
echo "pipeline done"
