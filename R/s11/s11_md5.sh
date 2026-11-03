#!/bin/bash
# s11_md5.sh -- S11 gate GC1: md5 of results_full/ and results_v4/ (excluding S11 and the S11 variants). Usage: s11_md5.sh before|after
set -u; cd /home/user/MMM; W=${1:?before or after}; O=results_v4/S11/md5_inputs_$W.txt; mkdir -p results_v4/S11
find results_full results_v4 -type f ! -path 'results_v4/S11/*' ! -path 'results_v4/variant_cv/*' ! -path 'results_v4/variant_t4cv/*' ! -path 'results_v4/variant_t2cv/*' -print0 | sort -z | xargs -0 -P 8 -n 200 md5sum | sort -k2 > $O
wc -l < $O
