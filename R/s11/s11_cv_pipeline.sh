#!/bin/bash
# s11_cv_pipeline.sh -- S11 (D-90, D-91): the S4-S6 chain of the table variant on the row set "cv" (variant nc with the
# D-90 no-result datasets not usable for the MMM with the convergence check). Same scripts as R/s9/s9_nc_pipeline.sh;
# the keys and the build step differ.
set -u
cd ~/MMM
export MMM_VARIANT=cv OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
L=reports/morris/variant_cv/logs; mkdir -p $L
for s in R/s11/s11_cv_keys.R R/s11/s11_cv_build.R R/postrun/agg_v4.R R/s5/s5_build.R R/s5/s5_run.R R/s5/s5_checks.R \
         R/s6/s6_run.R R/s6/s6_run_boby.R R/s6/s6_checks.R R/s6/s6_tables.R R/s6/s6_figdata.R R/s6/s6_dp_rates.R \
         R/s6/s6_rule.R R/s6/s6_rule_safety.R R/s6/s6_m1c_reasons.R R/s6/s6_m1c_glmer_extreme.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
echo "cv chain done"
