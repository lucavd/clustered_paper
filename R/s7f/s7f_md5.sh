#!/bin/bash
# s7f_md5.sh <before|after> -- GF1: md5 of results_full/ (against the S7b after-list) and of the row sets read in S7f
cd ~/MMM; tag=$1; O=results_v4/S7f; mkdir -p $O
( cd results_full && find . -name '*.rds' -type f | sort | xargs md5sum ) > $O/md5_results_full_$tag.txt
d=$(diff results_v4/S7b/md5_results_full_after.txt $O/md5_results_full_$tag.txt | grep -c '^[<>]' || true)
echo "tag=$tag results_full files=$(wc -l < $O/md5_results_full_$tag.txt) diff_vs_S7b_after=$d"
( for D in results_v4/S7d/rows results_v4/S7e/rows results_v4/S7d2/rows; do find $D -type f | sort | xargs md5sum; done; md5sum results_v4/variant_cc/S4/rows_v4.rds ) > $O/md5_inputs_$tag.txt
echo "tag=$tag input files=$(wc -l < $O/md5_inputs_$tag.txt) digest=$(md5sum < $O/md5_inputs_$tag.txt | cut -c1-32)"
