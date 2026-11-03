#!/bin/bash
# s6_g1.sh <tag> -- G1 of S6: md5 of results_full/, results_v4 (S2, S3, S4) and the S5 intermediates
set -u
TAG=$1; R=$HOME/MMM; LOG=$R/reports/morris/S6_logs/S6_g1_$TAG.log
mkdir -p $R/reports/morris/S6_logs
{
echo "G1 $TAG $(date '+%F %T')"
cd $R/results_full && n=$(md5sum -c --quiet ../baseline/results_full_md5.txt 2>&1 | wc -l); echo "G1 results_full: $(wc -l < ../baseline/results_full_md5.txt) checked, $n FAIL lines"
for s in S2 S3 S4; do cd $R/results_v4 && n=$(md5sum -c --quiet ../baseline/results_v4_${s}_md5.txt 2>&1 | wc -l); echo "G1 results_v4 $s: $(wc -l < ../baseline/results_v4_${s}_md5.txt) checked, $n FAIL lines"; done
cd $R && n=$(md5sum -c --quiet reports/morris/S5_logs/S5_intermediates_md5.txt 2>&1 | wc -l); echo "G1 results_v4 S5: $(wc -l < reports/morris/S5_logs/S5_intermediates_md5.txt) checked, $n FAIL lines"
} > $LOG
cat $LOG
