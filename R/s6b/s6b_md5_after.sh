#!/bin/bash
# S6b gates GA2/GB5: md5 of results_full/ and of the two S6cc inputs after the step, against the lists written before (s6b_thr_job.sh).
cd ~/MMM
( cd results_full && find . -name '*.rds' -type f | sort | xargs md5sum ) > results_v4/S6b/md5_results_full_after.txt
md5sum results_v4/S6cc/cc_all.rds results_v4/variant_cc/S4/rows_v4.rds > results_v4/S6b/md5_inputs_after.txt
n1=$(wc -l < results_v4/S6b/md5_results_full_before.txt); d1=$(diff results_v4/S6b/md5_results_full_before.txt results_v4/S6b/md5_results_full_after.txt | grep -c '^<' || true)
n2=$(wc -l < results_v4/S6b/md5_inputs_before.txt); d2=$(diff results_v4/S6b/md5_inputs_before.txt results_v4/S6b/md5_inputs_after.txt | grep -c '^<' || true)
printf "group,n,n_diff,same\nresults_full,%s,%s,%s\nS6cc inputs (cc_all.rds; variant_cc rows_v4.rds),%s,%s,%s\n" "$n1" "$d1" $([ "$d1" = 0 ] && echo TRUE || echo FALSE) "$n2" "$d2" $([ "$d2" = 0 ] && echo TRUE || echo FALSE) > tables/V4_S6b_md5_check.csv
cat tables/V4_S6b_md5_check.csv
