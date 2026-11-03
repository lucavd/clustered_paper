# REV2M_i.R -- which critical-value convention reproduces P_power_sizeadj.csv?
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); LRT <- "glmer_lrt"
load1 <- function(pref) {
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  one <- function(f) {
    x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
    x <- x[term=="time:trt"]; if (!nrow(x)) return(NULL)
    x[, model := as.integer(model)]
    x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
    x[, usable_strict := ifelse(method==LRT, converged & is.finite(p),
        converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
        abs(estimate)<10 & se>0 & se<10)]
    x[, usable_loose := converged & is.finite(p)]   # p-only filter
    x[usable_loose==TRUE, .(method,model,N,ntime,label,p,usable_strict)]
  }
  L <- mclapply(fs, one, mc.cores=45)
  rbindlist(L[!sapply(L,is.null)], fill=TRUE)
}
pub <- fread(path.expand("~/MMM/tables/P_power_sizeadj.csv"))
LB <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
for (pref in names(LB)) {
  lab <- LB[[pref]]; d <- load1(pref)
  cat("\n=========", lab, "=========\n")
  for (flt in c("strict","loose")) {
    dd <- if (flt=="strict") d[usable_strict==TRUE] else d
    for (lev in list(c("method","model","N","ntime"), c("method","N"), c("method","N","ntime"))) {
      crit <- dd[label=="manual", .(cv=quantile(p,0.05,type=1,names=FALSE)), by=lev]
      pw <- merge(dd[label=="mmlb"], crit, by=lev)
      agg <- pw[, .(n=.N, raw=100*mean(p<0.05), adj=100*mean(p<cv)), by=method]
      pb <- pub[dgp==lab, .(n_pub=sum(n), raw_pub=sum(power_raw*n)/sum(n), adj_pub=sum(power_adj*n)/sum(n)), by=method]
      mm <- merge(agg, pb, by="method")
      cat(sprintf("filter=%-6s cv-level=%-28s  max|adj-pub|=%.2f  max|raw-pub|=%.2f  max|n-npub|=%d\n",
        flt, paste(lev, collapse="+"), max(abs(mm$adj-mm$adj_pub)), max(abs(mm$raw-mm$raw_pub)),
        max(abs(mm$n-mm$n_pub))))
    }
  }
  # the best-matching convention printed in full
  dd <- d
  crit <- dd[label=="manual", .(cv=quantile(p,0.05,type=1,names=FALSE)), by=c("method","N")]
  pw <- merge(dd[label=="mmlb"], crit, by=c("method","N"))
  agg <- pw[, .(n=.N, raw=round(100*mean(p<0.05),2), adj=round(100*mean(p<cv),2)), by=method]
  pb <- pub[dgp==lab, .(n_pub=sum(n), raw_pub=round(sum(power_raw*n)/sum(n),2),
                        adj_pub=round(sum(power_adj*n)/sum(n),2)), by=method]
  print(as.data.frame(merge(agg, pb, by="method")[order(-adj)]), row.names=FALSE)
}
cat("\n########## 3.6 KC/uncorrected under both weightings\n")
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
S <- R[set=="own" & label=="manual"]
w <- S[, .(se_w=sum(mean_se*n)/sum(n), se_u=mean(mean_se)), by=.(method,N,dgp)]
for (dg in unique(w$dgp)) {
  x <- dcast(w[dgp==dg], N~method, value.var=c("se_w","se_u"))
  cat(sprintf("%s: KC/GEE replicate-weighted %.4f (N=10) -> %.4f (N=145); cell-unweighted %.4f -> %.4f\n",
    dg, x[N==10]$se_w_gee_kc/x[N==10]$se_w_gee, x[N==145]$se_w_gee_kc/x[N==145]$se_w_gee,
    x[N==10]$se_u_gee_kc/x[N==10]$se_u_gee, x[N==145]$se_u_gee_kc/x[N==145]$se_u_gee))
}
cat("\nH1 arm instead of H0:\n")
S1 <- R[set=="own" & label=="mmlb"]
w1 <- S1[, .(se_w=sum(mean_se*n)/sum(n), se_u=mean(mean_se)), by=.(method,N,dgp)]
for (dg in unique(w1$dgp)) {
  x <- dcast(w1[dgp==dg], N~method, value.var=c("se_w","se_u"))
  cat(sprintf("%s: KC/GEE weighted %.4f -> %.4f ; unweighted %.4f -> %.4f\n", dg,
    x[N==10]$se_w_gee_kc/x[N==10]$se_w_gee, x[N==145]$se_w_gee_kc/x[N==145]$se_w_gee,
    x[N==10]$se_u_gee_kc/x[N==10]$se_u_gee, x[N==145]$se_u_gee_kc/x[N==145]$se_u_gee))
}
cat("\nDONE_I\n")
