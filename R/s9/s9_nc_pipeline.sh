#!/bin/bash
# s9_nc_pipeline.sh -- S9 (D-82, D-83): the S4-S6 chain of the table variant on the row set "nc" (variant kcfc with the
# constant-outcome datasets not usable). Same scripts as R/s7f/s7f_kcfc_pipeline.sh; only the build step differs.
set -u
cd ~/MMM
export MMM_VARIANT=nc OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
L=reports/morris/variant_nc/logs; mkdir -p $L
for s in R/s9/s9_nc_build.R R/postrun/agg_v4.R R/s5/s5_build.R R/s5/s5_run.R R/s5/s5_checks.R \
         R/s6/s6_run.R R/s6/s6_run_boby.R R/s6/s6_checks.R R/s6/s6_tables.R R/s6/s6_figdata.R R/s6/s6_dp_rates.R \
         R/s6/s6_rule.R R/s6/s6_rule_safety.R R/s6/s6_m1c_reasons.R R/s6/s6_m1c_glmer_extreme.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
echo "nc chain done"
