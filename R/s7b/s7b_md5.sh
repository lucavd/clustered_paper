#!/bin/bash
# s7b_md5.sh <before|after> -- GC1: md5 of results_full/ (against results_v4/S7/md5_after.txt) and of the S7b inputs.
cd ~/MMM; tag=$1; O=results_v4/S7b; mkdir -p $O
( cd results_full && find . -name '*.rds' -type f | sort | xargs md5sum ) > $O/md5_results_full_$tag.txt
md5sum results_v4/variant_cc/S4/rows_v4.rds results_v4/S6b/ns/rows_ns.rds tables/variant_cc/V4_S6_rule_safety_dp.csv tables/V4_S6b_ns_rule_safety_dp.csv \
  tables/V4_S6b_ns_dp.csv tables/V4_S2_degenerate_designs.csv tables/V4_S4_truth_by_cell.csv tables/V4_S6b_ns_truth.csv rev2/T3_draws.rds rev2/T3_subsample_rerun.csv > $O/md5_inputs_$tag.txt
d0=$(diff <(sed 's/: OK$//' results_v4/S7/md5_after.txt | sort) <(awk '{print $2}' $O/md5_results_full_$tag.txt | sort) | grep -c '^[<>]' || true)
# results_v4/S7/md5_after.txt is a `md5sum -c` log (file: OK); compare file lists and then check sums against S6b's list
d1=$(diff results_v4/S6b/md5_results_full_after.txt $O/md5_results_full_$tag.txt | grep -c '^[<>]' || true)
echo "tag=$tag files=$(wc -l < $O/md5_results_full_$tag.txt) list_diff_vs_S7=$d0 sum_diff_vs_S6b=$d1"
