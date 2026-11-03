
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
fs <- list.files(RES, pattern="[.]rds$", full.names=TRUE)
cat("files:", length(fs), "\n")
KEEP <- c("gee","gee_kc","gee_md","gee_fg")
rd <- function(x) {
  z <- as.data.table(readRDS(x))
  z <- z[term=="time:trt" & method %in% KEEP,
         .(estimate, se, method, N, ntime, label, model, dgp, sim_id, converged)]
  z[, cid := sub("^.*_c([0-9]+)[.]rds$","\\1", basename(x))]
  z
}
D <- rbindlist(mclapply(fs, rd, mc.cores=40))
D[, model := as.integer(model)]
cat("rows read:", nrow(D), "\n")
W <- dcast(D, dgp+model+N+ntime+label+cid+sim_id ~ method,
           value.var="se", fun.aggregate=function(x) x[1])
cat("wide replicate rows:", nrow(W), "\n")
setnames(W, c("gee","gee_kc","gee_md","gee_fg"), c("se_gee","se_kc","se_md","se_fg"))

ratio_report <- function(num, nm) {
  S <- W[is.finite(se_gee) & se_gee > 0 & is.finite(get(num)) & get(num) > 0]
  S[, rr := get(num)/se_gee]
  qs <- quantile(S$rr, c(.01,.25,.5,.75,.95,.99,.999,.9999,1), type=7)
  cat("\n===== ", nm, " / uncorrected =====\n")
  cat("rows (both SE finite and > 0, NO usability filter):", nrow(S), "\n")
  print(signif(qs, 6))
  cat("max =", format(max(S$rr), digits=6), "  log10(max) =",
      round(log10(max(S$rr)),4), "\n")
  for (th in c(1e2,1e6,1e12,1e14,1e15,1e21)) {
    E <- S[rr > th]
    cat(sprintf("  > %-6g : rows %6d | cells %5d | design points %5d\n", th, nrow(E),
        uniqueN(E[, .(dgp,model,N,ntime,label)]), uniqueN(E[, .(dgp,model,N,ntime)])))
  }
  E <- S[rr > 1e12]
  if (nrow(E)) { cat("  ntime of rows > 1e12:\n"); print(table(E$ntime))
                 cat("  N of rows > 1e12:\n"); print(table(E$N))
                 cat("  dgp of rows > 1e12:\n"); print(table(E$dgp)) }
  E14 <- S[rr > 1e14]
  if (nrow(E14)) { cat("  rows > 1e14, ntime table:\n"); print(table(E14$ntime)) }
  invisible(S[, .(dgp,model,N,ntime,label,cid,sim_id,rr)][order(-rr)][1:20])
}
top_kc <- ratio_report("se_kc","Kauermann-Carroll")
top_md <- ratio_report("se_md","Mancl-DeRouen")
top_fg <- ratio_report("se_fg","Fay-Graubard")
cat("\n--- top 20 KC ratios ---\n"); print(top_kc)
cat("\n--- top 20 MD ratios ---\n"); print(top_md)
cat("\n--- top 20 FG ratios ---\n"); print(top_fg)

## does a usability filter change the picture? (declare both denominators)
U <- W[is.finite(se_gee) & se_gee>0 & se_gee<10 & is.finite(se_kc) & se_kc>0 & se_kc<10]
cat("\nrows surviving the usability bound (SE < 10) on both gee and gee_kc:", nrow(U),
    "  max ratio there =", format(max(U$se_kc/U$se_gee), digits=6), "\n")

out <- rbindlist(list(
  cbind(corr="kc", top_kc), cbind(corr="md", top_md), cbind(corr="fg", top_fg)))
fwrite(out, file.path(OUT,"R2_kc_ratio_top.csv"))
S <- W[is.finite(se_gee) & se_gee>0 & is.finite(se_kc) & se_kc>0]
S[, rr := se_kc/se_gee]
fwrite(S[, .(q = c("p01","p25","p50","p75","p95","p99","p999","p9999","max"),
             v = c(quantile(rr, c(.01,.25,.5,.75,.95,.99,.999,.9999), type=7), max(rr)))],
       file.path(OUT,"R2_kc_ratio_quantiles.csv"))
cat("\nDONE_C\n")
