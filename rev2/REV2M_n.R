# REV2M_n.R -- residual numbers: glmer empirical limit at N=145,T=10; which target
# agg_all5.R actually scored glmer against; §3.15 recomputation against plim.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method=="glmer" & label=="mmlb"]; if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10,
    .(dgp,model,N,ntime,estimate,se)]
}
D <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)), fill=TRUE)
cat("glmer usable H1 rows:", nrow(D), "\n")
cat("\n########## §3.15: 'mean estimate at 145 clusters and ten occasions is 0.279'\n")
z <- D[N==145 & ntime==10, .(n=.N, mean_est=round(mean(estimate),4), mcse=round(sd(estimate)/sqrt(.N),4)),
       by=.(dgp,model)][order(dgp,model)]
print(as.data.frame(z), row.names=FALSE)
cat("\npooled over the four main scenarios (M1,M2,M3,M5; M4 has a smaller effect):\n")
print(as.data.frame(D[N==145 & ntime==10 & model!=4,
  .(n=.N, mean_est=round(mean(estimate),4), mcse=round(sd(estimate)/sqrt(.N),4)), by=dgp]), row.names=FALSE)
P <- fread(path.expand("~/MMM/tables/O_glmer_plim.csv"))
cat("\nO_glmer_plim.csv at ntime=10:\n")
print(as.data.frame(P[ntime==10][, .(dgp,model,plim=round(plim,4),n_plim)]), row.names=FALSE)
cat("\n########## which target did agg_all5.R score glmer against?\n")
cat("agg_all5.R lines 89-90: t_gen <- t_marg ; t_plim <- ifelse(method=='glmer' & is.finite(plim), plim, t_marg)\n")
cat("=> cov_gen for glmer uses the MARGINAL target (0.200), not the conditional 0.232.\n")
O <- fread(path.expand("~/MMM/tables/O_coverage_all.csv"))
og <- O[set=="own" & method=="glmer", .(n=sum(n), cov_gen=round(sum(cov_gen*n)/sum(n),2),
                                        cov_plim=round(sum(cov_plim*n)/sum(n),2)), by=dgp]
print(as.data.frame(og), row.names=FALSE)
cat("\n########## coverage recomputed against the plim, all methods\n")
zq <- qnorm(0.975)
PL <- P[, .(dgp,model,ntime,plim)]
d2 <- merge(D, PL, by=c("dgp","model","ntime"), all.x=TRUE)
d2[, cov_plim := abs(estimate-plim)<=zq*se]
d2[, cov_cond := abs(estimate-ifelse(model==4,0.058,0.232))<=zq*se]
d2[, cov_marg := abs(estimate-ifelse(model==4,0.050,0.200))<=zq*se]
print(as.data.frame(d2[, .(n=.N, vs_plim=round(100*mean(cov_plim),2), vs_cond=round(100*mean(cov_cond),2),
   vs_marg=round(100*mean(cov_marg),2)), by=dgp]), row.names=FALSE)
cat("\n########## §2.4 filter wording: does the code use 0 < SE < 10 or SE < 10?\n")
print(grep("se >|se <|abs\\(estimate\\)", readLines(path.expand("~/MMM/R/postrun/agg_all5.R")), value=TRUE))
cat("\nDONE_N\n")
