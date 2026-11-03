#!/bin/bash
# Installs fewclust from commit c688ecc (branch package) into ~/R/s9_fewclust_lib and runs the citation condition.
set -euo pipefail
C=c688eccdee52026b79370b107bf3a2e852f8d08c; LIB=$HOME/R/s9_fewclust_lib; T=$(mktemp -d)
mkdir -p $LIB
git -C /home/user/MMM-package archive $C fewclust | tar -x -C $T
R_LIBS=$LIB:/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu R CMD INSTALL -l $LIB $T/fewclust > $LIB/install.log 2>&1
echo $C > $LIB/fewclust_commit.txt
cd /home/user/MMM && Rscript R/s9/s9_fewclust_condition.R
rm -rf $T
