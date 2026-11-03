# REV2M_j.R -- the KC SE ratio computed PAIRED per replicate (as kc.log does),
# the extreme-ratio claim ("order 10^21"), and the 1.344 -> 0.813 sequence.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method %in% c("gee","gee_kc","gee_md","gee_fg")]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  w <- dcast(x[converged==TRUE & is.finite(se) & se>0, .(dgp,model,N,ntime,label,sim_id,method,se,estimate)],
             dgp+model+N+ntime+label+sim_id~method, value.var=c("se","estimate"))
  setnames(w, sub("^se_","",names(w)))
  if (!all(c("gee","gee_kc") %in% names(w))) return(NULL)
  w
}
fs <- list.files(RES, full.names=TRUE)
L <- mclapply(fs, one, mc.cores=45)
W <- rbindlist(L[!sapply(L,is.null)], fill=TRUE)
cat("paired replicate rows with both gee and gee_kc SE:", nrow(W[!is.na(gee) & !is.na(gee_kc)]), "\n")
W <- W[!is.na(gee) & gee>0]
W[, `:=`(r_kc=gee_kc/gee, r_md=gee_md/gee, r_fg=gee_fg/gee)]
cat("\n===== paired ratio to the uncorrected sandwich, NO |estimate|/se<10 filter =====\n")
cat(sprintf("KC: median %.4f  mean %.4f  max %.4g  min %.4g  (n=%d)\n",
  median(W$r_kc,na.rm=TRUE), mean(W$r_kc,na.rm=TRUE), max(W$r_kc,na.rm=TRUE), min(W$r_kc,na.rm=TRUE),
  sum(is.finite(W$r_kc))))
cat(sprintf("MD: median %.4f  max %.4g | FG: median %.4f  max %.4g\n",
  median(W$r_md,na.rm=TRUE), max(W$r_md,na.rm=TRUE), median(W$r_fg,na.rm=TRUE), max(W$r_fg,na.rm=TRUE)))
cat("\nreplicates with KC/uncorrected above 100 / 1e6 / 1e12 / 1e20:\n")
for (thr in c(1e2,1e6,1e12,1e20)) cat(sprintf("  > %.0e : %d\n", thr, sum(W$r_kc>thr, na.rm=TRUE)))
cat(sprintf("largest 5 KC ratios: %s\n", paste(round(sort(W$r_kc, decreasing=TRUE)[1:5],4), collapse=", ")))
cat("\n===== paired mean ratio by cluster count (the 1.344 -> 0.813 claim) =====\n")
byN <- W[, .(n=.N, kc=round(mean(r_kc,na.rm=TRUE),4), kc_med=round(median(r_kc,na.rm=TRUE),4),
             md=round(mean(r_md,na.rm=TRUE),4), fg=round(mean(r_fg,na.rm=TRUE),4)), by=N][order(N)]
print(as.data.frame(byN), row.names=FALSE)
cat("\n===== GenBinaryY only, MEDIAN paired ratio by N (unfiltered) =====\n")
print(as.data.frame(W[dgp=="GenBinaryY", .(n=.N, kc_med=round(median(r_kc,na.rm=TRUE),4),
  md_med=round(median(r_md,na.rm=TRUE),4), fg_med=round(median(r_fg,na.rm=TRUE),4)), by=N][order(N)]), row.names=FALSE)
cat("\n===== filtered set (se<10, |est|<10) =====\n")
Wf <- W[gee<10 & gee_kc<10 & abs(estimate_gee)<10]
byN2 <- Wf[, .(n=.N, kc=round(mean(r_kc,na.rm=TRUE),4), md=round(mean(r_md,na.rm=TRUE),4),
               fg=round(mean(r_fg,na.rm=TRUE),4)), by=N][order(N)]
print(as.data.frame(byN2), row.names=FALSE)
cat(sprintf("filtered paired MEAN: KC %.4f at N=10 -> %.4f at N=145 ; overall %.4f\n",
  byN2[N==10]$kc, byN2[N==145]$kc, mean(Wf$r_kc, na.rm=TRUE)))
byN2m <- Wf[, .(n=.N, kc=round(median(r_kc,na.rm=TRUE),4), md=round(median(r_md,na.rm=TRUE),4),
                fg=round(median(r_fg,na.rm=TRUE),4)), by=N][order(N)]
cat("filtered MEDIAN by N:\n"); print(as.data.frame(byN2m), row.names=FALSE)
cat(sprintf("filtered median: KC %.4f -> %.4f ; MD %.4f -> %.4f ; FG %.4f -> %.4f\n",
  byN2m[N==10]$kc, byN2m[N==145]$kc, byN2m[N==10]$md, byN2m[N==145]$md,
  byN2m[N==10]$fg, byN2m[N==145]$fg))
cat(sprintf("GenBinaryY filtered median KC: %.4f -> %.4f\n",
  median(Wf[dgp=="GenBinaryY" & N==10]$r_kc, na.rm=TRUE), median(Wf[dgp=="GenBinaryY" & N==145]$r_kc, na.rm=TRUE)))
cat(sprintf("unfiltered paired: KC %.4f -> %.4f ; overall %.3f\n",
  byN[N==10]$kc, byN[N==145]$kc, mean(W$r_kc, na.rm=TRUE)))
cat("\n===== ratio of MEAN SEs (not mean of ratios), filtered, by N =====\n")
byN3 <- Wf[, .(n=.N, kc=round(mean(gee_kc)/mean(gee),4), md=round(mean(gee_md,na.rm=TRUE)/mean(gee),4),
               fg=round(mean(gee_fg,na.rm=TRUE)/mean(gee),4)), by=N][order(N)]
print(as.data.frame(byN3), row.names=FALSE)
cat(sprintf("ratio-of-means: KC %.4f at N=10 -> %.4f at N=145\n", byN3[N==10]$kc, byN3[N==145]$kc))
cat("\nDONE_J\n")
