#!/bin/bash
# s7g_md5.sh <before|after> -- gate 1: md5 of results_full/ against the S7b after-list (as R/s7e/s7e_md5.sh)
cd ~/MMM; tag=$1; O=results_v4/S7g; mkdir -p $O
( cd results_full && find . -name '*.rds' -type f | sort | xargs md5sum ) > $O/md5_results_full_$tag.txt
d=$(diff results_v4/S7b/md5_results_full_after.txt $O/md5_results_full_$tag.txt | grep -c '^[<>]' || true)
echo "tag=$tag files=$(wc -l < $O/md5_results_full_$tag.txt) diff_vs_S7b_after=$d"
