
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
F1 <- fread(path.expand("~/MMM/tables/F_truth_ar1.csv"))
R2 <- fread(path.expand("~/MMM/rev2/R2_ar1_cond_true.csv"))
cat("F_truth_ar1: righe", nrow(F1), " N distinti:", paste(sort(unique(F1$N)),collapse=","), "\n")
cat("R2_ar1_cond_true: righe", nrow(R2), " indicizzato per (model,ntime) solo\n\n")
# la media di F1 su N per (model,ntime) coincide con R2?
a <- F1[, .(f_mean = mean(cond_used)), by=.(model,ntime)]
b <- merge(a, R2[, .(model, ntime, r = t_cond_true)], by=c("model","ntime"))
cat("massimo scarto fra la media di F_truth_ar1 su N e R2_ar1_cond_true:",
    format(max(abs(b$f_mean - b$r)), digits=4), "\n")
cat("scarto mediano:", format(median(abs(b$f_mean - b$r)), digits=4), "\n")
print(head(b[order(-abs(f_mean-r))], 5))
