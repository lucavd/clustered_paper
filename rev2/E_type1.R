
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
fs <- list.files(RES, pattern="_manual_c[0-9]+[.]rds$", full.names=TRUE)
k <- sub("_c[0-9]+[.]rds$","",basename(fs)); fs <- fs[k %in% names(which(table(k)==10))]
cat("H0 files:", length(fs), "\n")
D <- rbindlist(mclapply(fs, function(x){ z <- as.data.table(readRDS(x))
  z[term=="time:trt", .(estimate,se,p,method,N,ntime,model,converged,dgp)] }, mc.cores=40))
D[, model := as.integer(model)]
H0 <- rbind(D[method %in% WALD & converged & is.finite(estimate) & is.finite(se) &
              abs(estimate)<10 & se>0 & se<10 & is.finite(p)],
            D[method=="glmer_lrt" & converged & is.finite(p)])
cat("\n### TYPE I ERROR at nominal 5%, pooled over usable H0 replicates (row-weighted)\n")
print(H0[, .(n=.N, typeI=round(100*mean(p<0.05),4),
             mcse=round(100*sqrt(mean(p<0.05)*(1-mean(p<0.05))/.N),4)),
         by=.(dgp,method)][order(dgp,-typeI)])
fwrite(H0[, .(n=.N, typeI=100*mean(p<0.05)), by=.(dgp,method)],
       file.path(OUT,"R2_typeI_pooled.csv"))
cat("\nDONE_E\n")
