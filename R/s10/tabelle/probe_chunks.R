
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
x <- readRDS("/home/user/MMM/results_v4/S7d/rows/GenBinaryY_m1c_N10_T2_manual_c01.rds")
str(x$recs[1:2], vec.len=2); str(x$rows[1:3], vec.len=2)
cat("---- S7e\n")
y <- readRDS("/home/user/MMM/results_v4/S7e/rows/copula_m1_N105_T10_manual_c01.rds")
print(names(y)); for (n in names(y)) {cat("==",n,"\n"); if (is.data.frame(y[[n]])) str(as.data.frame(y[[n]])[1:2,], vec.len=2) else str(y[[n]], max.level=1)}
