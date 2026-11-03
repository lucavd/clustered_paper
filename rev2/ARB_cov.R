
# Arbitra 93,16 (track Metodi) contro 93,19 (track Tabelle): copertura del GLMM su
# GLMM-AR1 contro la verita' condizionale calibrata per (scenario, T).
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
TR  <- fread(path.expand("~/MMM/rev2/R2_ar1_cond_true.csv"))
setnames(TR, "t_cond_true", "tcond")
fs <- list.files(RES, pattern="^GLMMAR1_.*_mmlb_c[0-9]+[.]rds$")
cat("chunk H1 AR1:", length(fs), "\n")
zq <- qnorm(0.975)
one <- function(f) {
  x <- try(as.data.table(readRDS(file.path(RES,f))), silent=TRUE)
  if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method=="glmer"]
  if (!nrow(x)) return(NULL)
  x[, converged := as.logical(converged)]; x[is.na(converged), converged := FALSE]
  u <- x[converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
         abs(estimate) < 10 & se > 0 & se < 10]
  if (!nrow(u)) return(NULL)
  u[, .(model=as.integer(model), ntime, estimate, se)]
}
L <- mclapply(fs, one, mc.cores=40)
D <- rbindlist(Filter(is.data.table, L))
cat("repliche utilizzabili:", nrow(D), "\n")
D <- merge(D, TR[, .(model, ntime, tcond)], by=c("model","ntime"), all.x=TRUE)
cat("senza verita' condizionale:", D[is.na(tcond), .N], "\n")
D[, cov_cond := abs(estimate - tcond) <= zq*se]
D[, cov_marg := abs(estimate - 0.200) <= zq*se]
cat(sprintf("\nCOPERTURA contro il condizionale per (scenario,T): %.4f  su n=%d\n",
            100*mean(D$cov_cond), nrow(D)))
cat(sprintf("COPERTURA contro il marginale 0.200:               %.4f\n", 100*mean(D$cov_marg)))
# scenario 4 ha bersaglio marginale 0.05: rifare con il marginale corretto per scenario
D[, tmarg := ifelse(model==4, 0.050, 0.200)]
D[, cov_marg2 := abs(estimate - tmarg) <= zq*se]
cat(sprintf("COPERTURA contro il marginale per scenario:        %.4f\n", 100*mean(D$cov_marg2)))
cat("\nper scenario (condizionale):\n")
print(D[, .(n=.N, cov=round(100*mean(cov_cond),4)), by=model][order(model)])
