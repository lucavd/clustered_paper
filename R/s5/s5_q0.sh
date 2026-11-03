#!/bin/bash
# s5_q0.sh <tag> -- G1 of S5: md5 of results_full/, results_v4 (S2, S3, S4 intermediates)
set -u
TAG=$1; R=$HOME/MMM; LOG=$R/reports/morris/S5_logs/S5_g1_$TAG.log
mkdir -p $R/reports/morris/S5_logs
{
echo "G1 $TAG $(date '+%F %T')"
cd $R/results_full && n=$(md5sum -c --quiet ../baseline/results_full_md5.txt 2>&1 | wc -l); echo "G1 results_full: $(wc -l < ../baseline/results_full_md5.txt) checked, $n FAIL lines"
cd $R/results_v4 && n=$(md5sum -c --quiet ../baseline/results_v4_S2_md5.txt 2>&1 | wc -l); echo "G1 results_v4 S2: $(wc -l < ../baseline/results_v4_S2_md5.txt) checked, $n FAIL lines"
cd $R/results_v4 && n=$(md5sum -c --quiet ../baseline/results_v4_S3_md5.txt 2>&1 | wc -l); echo "G1 results_v4 S3: $(wc -l < ../baseline/results_v4_S3_md5.txt) checked, $n FAIL lines"
cd $R/results_v4 && n=$(md5sum -c --quiet ../baseline/results_v4_S4_md5.txt 2>&1 | wc -l); echo "G1 results_v4 S4: $(wc -l < ../baseline/results_v4_S4_md5.txt) checked, $n FAIL lines"
} > $LOG
cat $LOG
