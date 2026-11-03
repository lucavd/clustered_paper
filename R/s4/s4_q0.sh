#!/bin/bash
# s4_q0.sh <tag> -- Q0 of S4: md5 of results_full/, results_v4 (S2, S3) and rows_timetrt.rds
set -u
TAG=$1; R=$HOME/MMM; LOG=$R/reports/morris/S4_logs/S4_q0_$TAG.log
mkdir -p $R/reports/morris/S4_logs
{
echo "Q0 $TAG $(date '+%F %T')"
cd $R/results_full && n=$(md5sum -c --quiet ../baseline/results_full_md5.txt 2>&1 | wc -l); echo "Q0 results_full: $(wc -l < ../baseline/results_full_md5.txt) checked, $n FAIL lines"
cd $R/results_v4 && n=$(md5sum -c --quiet ../baseline/results_v4_S2_md5.txt 2>&1 | wc -l); echo "Q0 results_v4 S2: $(wc -l < ../baseline/results_v4_S2_md5.txt) checked, $n FAIL lines"
cd $R/results_v4 && n=$(md5sum -c --quiet ../baseline/results_v4_S3_md5.txt 2>&1 | wc -l); echo "Q0 results_v4 S3: $(wc -l < ../baseline/results_v4_S3_md5.txt) checked, $n FAIL lines"
echo "Q0 rows_timetrt.rds $(md5sum $R/results_v4/S2/rows_timetrt.rds | cut -d' ' -f1) (expected 092a5ad60b60e0ea28cb6b78ba2adbeb)"
} > $LOG
cat $LOG
