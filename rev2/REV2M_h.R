# REV2M_h.R -- F (size-adjusted power reproduced from the chunks), 3.6 SE ratios,
#              provenance of the hand-computation logs, run-log error accounting.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink"); LRT <- "glmer_lrt"

cat("########## G-provenance: do kca.log / fga.log / mda.log exist anywhere?\n")
print(system("find ~ -name 'kca.log' -o -name 'fga.log' -o -name 'mda.log' -o -name 'kc*.log' -o -name 'fg*.log' 2>/dev/null | head -20", intern=TRUE))
cat("R scripts that CALL geesmv outside the pipeline:\n")
print(system("grep -rln 'GEE.var' ~/MMM --include='*.R' 2>/dev/null | grep -v renv", intern=TRUE))
cat("any script computing a leverage/hat matrix by hand:\n")
print(system("grep -rln 'solve(diag\\|diag(n_i)\\|I - H\\|leverage\\|msqrt\\|sqrtm' ~/MMM --include='*.R' 2>/dev/null | grep -v renv", intern=TRUE))

cat("\n########## 3.6: KC SE / uncorrected SE by cluster count, from REV_cells.csv\n")
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
S <- R[set=="own" & label=="manual", .(mean_se=sum(mean_se*n)/sum(n), n=sum(n)), by=.(method,N,dgp)]
W <- dcast(S, N+dgp~method, value.var="mean_se")
W[, `:=`(kc_over_gee=round(gee_kc/gee,4), md_over_gee=round(gee_md/gee,4), fg_over_gee=round(gee_fg/gee,4))]
print(as.data.frame(W[dgp=="GenBinaryY", .(N,kc_over_gee,md_over_gee,fg_over_gee)]), row.names=FALSE)
cat("pooled over mechanisms:\n")
P <- S[, .(mean_se=sum(mean_se*n)/sum(n)), by=.(method,N)]
WP <- dcast(P, N~method, value.var="mean_se")
WP[, `:=`(kc=round(gee_kc/gee,4), md=round(gee_md/gee,4), fg=round(gee_fg/gee,4))]
print(as.data.frame(WP[, .(N,kc,md,fg)]), row.names=FALSE)
cat(sprintf("KC/uncorrected: %.4f at N=10, %.4f at N=145 (pooled); GenBinaryY only: %.4f -> %.4f\n",
  WP[N==10]$kc, WP[N==145]$kc, W[N==10 & dgp=="GenBinaryY"]$kc_over_gee,
  W[N==145 & dgp=="GenBinaryY"]$kc_over_gee))

cat("\n########## F: size-adjusted power recomputed per CELL from the chunks\n")
do_dgp <- function(pref) {
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  one <- function(f) {
    x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
    x <- x[term=="time:trt"]; if (!nrow(x)) return(NULL)
    x[, model := as.integer(model)]
    x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
    x[, usable := ifelse(method==LRT, converged & is.finite(p),
        converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
        abs(estimate)<10 & se>0 & se<10)]
    x[usable==TRUE, .(method,model,N,ntime,label,p)]
  }
  L <- mclapply(fs, one, mc.cores=45)
  d <- rbindlist(L[!sapply(L,is.null)], fill=TRUE)
  # per cell: critical value = 5th percentile of the H0 p-values; power = share of
  # H1 p-values below it
  crit <- d[label=="manual", .(cv = quantile(p, 0.05, type=1, names=FALSE), n0=.N),
            by=.(method,model,N,ntime)]
  pw <- merge(d[label=="mmlb"], crit, by=c("method","model","N","ntime"))
  cell <- pw[, .(n1=.N, raw=100*mean(p<0.05), adj=100*mean(p<cv)), by=.(method,model,N,ntime)]
  list(byN = cell[, .(n=sum(n1), raw=round(sum(raw*n1)/sum(n1),2), adj=round(sum(adj*n1)/sum(n1),2)),
                  by=.(method,N)],
       pooled = cell[, .(n=sum(n1), raw=round(sum(raw*n1)/sum(n1),2),
                         adj=round(sum(adj*n1)/sum(n1),2)), by=method][order(-adj)])
}
pub <- fread(path.expand("~/MMM/tables/P_power_sizeadj.csv"))
for (pref in c("GenBinaryY","GLMMnative","GLMMAR1")) {
  lab <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")[[pref]]
  cat("\n=========", lab, "=========\n")
  o <- do_dgp(pref)
  cat("-- recomputed, pooled over the grid (cell-level critical value) --\n")
  print(as.data.frame(o$pooled), row.names=FALSE)
  cat("-- published table, same pooling --\n")
  pp <- pub[dgp==lab, .(n=sum(n), raw=round(sum(power_raw*n)/sum(n),2),
                        adj=round(sum(power_adj*n)/sum(n),2)), by=method][order(-adj)]
  print(as.data.frame(pp), row.names=FALSE)
  m <- merge(o$pooled, pp, by="method", suffixes=c("_rec","_pub"))
  m[, `:=`(d_adj=round(adj_rec-adj_pub,2), d_raw=round(raw_rec-raw_pub,2), d_n=n_rec-n_pub)]
  cat("-- difference (recomputed minus published) --\n")
  print(as.data.frame(m[, .(method,adj_rec,adj_pub,d_adj,raw_rec,raw_pub,d_raw,d_n)]), row.names=FALSE)
  cat(sprintf("recomputed adjusted-power spread across the 8 methods: %.2f - %.2f\n",
              min(o$pooled$adj), max(o$pooled$adj)))
  saveRDS(o, path.expand(paste0("~/MMM/rev2/REV2M_power_", pref, ".rds")))
}

cat("\n########## run-log accounting for '~zero errors~'\n")
print(system("grep -c ERROR ~/MMM/rerun_m3m5.log ~/MMM/_superseded/monitoring/*.log 2>/dev/null", intern=TRUE))
print(system("grep -h 'DONE:\\|total chunks' ~/MMM/rerun_m3m5.log ~/MMM/_superseded/monitoring/*.log 2>/dev/null | head -20", intern=TRUE))
cat("\nDONE_H\n")
