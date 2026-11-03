#!/bin/bash
# s11_thr_resume.sh -- S11 (F13): rest of the S6 chain of a threshold variant after s6_run (s6_figdata was missing from
# R/s11/s11_thr_pipeline.sh; s6_dp_rates reads its nested-loop table). Usage: s11_thr_resume.sh t4cv|t2cv
set -u
cd ~/MMM
TAG=${1:?t4cv or t2cv}
export MMM_VARIANT=$TAG OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 S5_CORES=20 R_LIBS=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
L=reports/morris/variant_$TAG/logs; mkdir -p $L
for s in R/s6/s6_figdata.R R/s6/s6_dp_rates.R R/s6/s6_rule.R R/s6/s6_rule_safety.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
echo "$TAG resume done"
