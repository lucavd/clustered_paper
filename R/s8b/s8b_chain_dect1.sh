#!/bin/bash
# s8b_chain_dect1.sh -- DEC-T1 (D-73): rerun the kcfc S5-S6 chain from s5_build.R with the declared size-adjusted rule,
# then s7_extras.R and the independent Table 5 recomputation. Logs in reports/morris/S8b_logs/kcfc/.
set -u
cd ~/MMM
export MMM_VARIANT=kcfc OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 S5_CORES=${S5_CORES:-6}
L=reports/morris/S8b_logs/kcfc; mkdir -p $L
for s in R/s5/s5_build.R R/s5/s5_run.R R/s5/s5_checks.R \
         R/s6/s6_run.R R/s6/s6_run_boby.R R/s6/s6_checks.R R/s6/s6_tables.R R/s6/s6_figdata.R R/s6/s6_dp_rates.R \
         R/s6/s6_rule.R R/s6/s6_rule_safety.R R/s6/s6_m1c_reasons.R R/s6/s6_m1c_glmer_extreme.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
unset MMM_VARIANT
for s in R/s7/s7_extras.R R/s7f/s7f_table5_recompute.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
echo "chain done"
