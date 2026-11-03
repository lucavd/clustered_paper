#!/bin/bash
# s6cc_pipeline.sh -- D-32: the S4-S6 chain on the variant row set (MMM with convergence check in M3-M5).
# Same scripts as the study, with MMM_VARIANT=cc (paths redirected in R/s4/s4_common.R, R/s5/s5_common.R,
# R/s6/s6_common.R). Stops at the first failing step.
set -u
cd ~/MMM
export MMM_VARIANT=cc
L=reports/morris/variant_cc/logs; mkdir -p $L
for s in R/s6/s6cc_build.R R/postrun/agg_v4.R R/s5/s5_build.R R/s5/s5_run.R R/s5/s5_checks.R \
         R/s6/s6_run.R R/s6/s6_run_boby.R R/s6/s6_checks.R R/s6/s6_tables.R R/s6/s6_figdata.R R/s6/s6_dp_rates.R \
         R/s6/s6_rule.R R/s6/s6_rule_safety.R R/s6/s6_m1c_reasons.R R/s6/s6_m1c_glmer_extreme.R R/s6/s6_mmm_diag.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -5 $L/$b.log; exit 1; fi
done
echo "pipeline done"
