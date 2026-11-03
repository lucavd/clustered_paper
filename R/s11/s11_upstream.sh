#!/bin/bash
# s11_upstream.sh -- S11 (D-90, D-91): per-point aggregations upstream of the closing pipeline that read the MMM with the
# convergence check, in dependency order: M6-M7 tables -> S7b (t with N - p and N - 2; null shift) -> S7d -> S7e (fourth
# mechanism) and its summary; then the check-versus-study cells of Section 3.2 and the shares of interventions.
# The S7d2 scripts (sign-flip and its GEE comparison) do not read the marginalized model and are not rerun.
# Logs in reports/morris/S11_logs/upstream/.
set -u
cd ~/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
L=reports/morris/S11_logs/upstream; mkdir -p $L
for s in R/s9/s9_ns_tables.R R/s7b/s7b_df.R R/s7b/s7b_null.R R/s7d/s7d_agg.R R/s7e/s7e_agg.R R/s7e/s7e_summary.R \
         R/s11/s11_cc_cells.R R/s11/s11_cc_shares.R; do
  b=$(basename $s .R); t0=$(date +%s)
  if Rscript --vanilla $s > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi
done
echo "upstream done"
