# REV2M_o.R -- §3.13: glmer bias under both targets and both pooling conventions,
# and the "42.8% of the difference is the scoring convention" claim.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & label=="mmlb" & method %in% c("glmer","mmm","gee_fg")]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10,
    .(dgp,method,model,N,ntime,estimate,se)]
}
D <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)), fill=TRUE)
D[, t_cond := ifelse(model==4, 0.058, 0.232)]
D[, t_marg := ifelse(model==4, 0.050, 0.200)]
P <- fread(path.expand("~/MMM/tables/O_glmer_plim.csv"))[, .(dgp,model,ntime,plim)]
D <- merge(D, P, by=c("dgp","model","ntime"), all.x=TRUE)
G <- D[method=="glmer"]
cat("glmer rows:", nrow(G), "\n\n")
cat("===== mean SIGNED bias, replicate-weighted over the whole grid =====\n")
print(as.data.frame(G[, .(n=.N, vs_cond=round(mean(estimate-t_cond),4),
  vs_marg=round(mean(estimate-t_marg),4), vs_plim=round(mean(estimate-plim),4)), by=dgp]), row.names=FALSE)
cat(sprintf("pooled: cond %.4f | marg %.4f | plim %.4f\n",
  mean(G$estimate-G$t_cond), mean(G$estimate-G$t_marg), mean(G$estimate-G$plim, na.rm=TRUE)))
cat("\n===== mean of per-DESIGN-POINT bias (unweighted over 1440 points) =====\n")
dp <- G[, .(b_cond=mean(estimate-t_cond), b_marg=mean(estimate-t_marg),
            b_plim=mean(estimate-plim), n=.N), by=.(dgp,model,N,ntime)]
print(as.data.frame(dp[, .(points=.N, vs_cond=round(mean(b_cond),4), vs_marg=round(mean(b_marg),4),
   vs_plim=round(mean(b_plim,na.rm=TRUE),4)), by=dgp]), row.names=FALSE)
cat(sprintf("pooled design-point mean: cond %.4f | marg %.4f | plim %.4f\n",
  mean(dp$b_cond), mean(dp$b_marg), mean(dp$b_plim, na.rm=TRUE)))
cat("\n===== the '42.8% of the difference is the scoring convention' claim =====\n")
cat("interpretation A: (bias_marg - bias_cond)/bias_marg\n")
for (lab in c("replicate-weighted","design-point mean")) {
  bc <- if (lab=="replicate-weighted") mean(G$estimate-G$t_cond) else mean(dp$b_cond)
  bm <- if (lab=="replicate-weighted") mean(G$estimate-G$t_marg) else mean(dp$b_marg)
  cat(sprintf("  %-20s bias_cond %.4f  bias_marg %.4f  share %.1f%%\n", lab, bc, bm, 100*(bm-bc)/bm))
}
cat("interpretation B: relative to the comparator's bias (mmm)\n")
M <- D[method=="mmm"]; bm_mmm <- mean(M$estimate-M$t_marg)
cat(sprintf("  mmm bias vs marginal %.4f ; glmer-minus-mmm gap: marg %.4f, cond %.4f, share %.1f%%\n",
  bm_mmm, mean(G$estimate-G$t_marg)-bm_mmm, mean(G$estimate-G$t_cond)-bm_mmm,
  100*((mean(G$estimate-G$t_marg)-bm_mmm)-(mean(G$estimate-G$t_cond)-bm_mmm))/(mean(G$estimate-G$t_marg)-bm_mmm)))
cat("\n===== per-cell ABSOLUTE bias (the 31.9% claim) =====\n")
cat(sprintf("mean |bias| per design point: cond %.4f | marg %.4f | share %.1f%%\n",
  mean(abs(dp$b_cond)), mean(abs(dp$b_marg)), 100*(mean(abs(dp$b_marg))-mean(abs(dp$b_cond)))/mean(abs(dp$b_marg))))
cat("\n===== coverage gap (the '1.0 pp of a 4.1 pp gap' claim) =====\n")
zq <- qnorm(0.975)
cv <- D[, .(cov_cond=100*mean(abs(estimate-t_cond)<=zq*se),
            cov_marg=100*mean(abs(estimate-t_marg)<=zq*se), n=.N), by=.(method,dgp)]
print(as.data.frame(cv[, lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)]), row.names=FALSE)
gg <- cv[method=="glmer"]; ff <- cv[method=="gee_fg"]
cat(sprintf("glmer vs gee_fg coverage gap: on the marginal target %.2f pp | on glmer's conditional %.2f pp | difference %.2f pp\n",
  mean(ff$cov_marg)-mean(gg$cov_marg), mean(ff$cov_marg)-mean(gg$cov_cond),
  (mean(ff$cov_marg)-mean(gg$cov_marg))-(mean(ff$cov_marg)-mean(gg$cov_cond))))
cat("\nDONE_O\n")
