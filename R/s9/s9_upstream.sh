#!/bin/bash
# s9_upstream.sh -- S9 (D-82, D-83): per-point aggregations upstream of the closing pipeline, with the constant-outcome
# rule, in dependency order: M6-M7 tables -> S7b (t with N - p and N - 2; null shift) -> S7d -> S7d2 (sign-flip interval
# and its GEE comparison) -> S7e (fourth mechanism) and its summary. Logs in reports/morris/S9_logs/upstream/.
set -u
cd ~/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
L=reports/morris/S9_logs/upstream; mkdir -p $L
for s in R/s9/s9_ns_tables.R R/s7b/s7b_df.R R/s7b/s7b_null.R R/s7d/s7d_agg.R R/s7d2/s7d2_agg.R R/s7d2/s7d2_explain.R \
         R/s7d2/s7d2_width_groups.R R/s7e/s7e_agg.R R/s7e/s7e_summary.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
echo "upstream done"
