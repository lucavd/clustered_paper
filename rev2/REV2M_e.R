# REV2M_e.R -- E (usability filter, per method), C (coverage with T-varying truth),
#              D (grid counts), F (size-adjusted power: does the code exist?)
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(data.table); library(dplyr); library(tidyr); library(fastGHQuad)})
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
R[, model := as.integer(model)]

cat("########## D. grid counts from REV_cells.csv (built from the raw chunks)\n")
cat(sprintf("design points (dgp,model,N,ntime): %d\n", uniqueN(R[,.(dgp,model,N,ntime)])))
cat(sprintf("cells (design point x hypothesis): %d\n", uniqueN(R[,.(dgp,model,N,ntime,label)])))
cat(sprintf("N levels: %d  T levels: %d  models: %d  dgps: %d\n",
  uniqueN(R$N), uniqueN(R$ntime), uniqueN(R$model), uniqueN(R$dgp)))
cat("distinct n_attempt per cell (replicates actually present):\n")
print(R[set=="own", .(cells=uniqueN(paste(dgp,model,N,ntime,label))), by=n_attempt][order(-cells)][1:10])
att <- unique(R[set=="own", .(dgp,model,N,ntime,label,n_attempt)])
cat(sprintf("cells with n_attempt == 1000: %d of %d (%.2f%%)\n",
  sum(att$n_attempt==1000), nrow(att), 100*mean(att$n_attempt==1000)))
cat(sprintf("n_attempt: min %d  median %g  mean %.1f  max %d\n",
  min(att$n_attempt), median(att$n_attempt), mean(att$n_attempt), max(att$n_attempt)))
cat("cells with n_attempt < 1000, by T:\n")
print(att[n_attempt<1000, .(cells=.N, mean_attempt=round(mean(n_attempt),1),
                            min_attempt=min(n_attempt)), by=ntime][order(ntime)])
cat("\nby dgp:\n"); print(att[, .(cells=.N, mean_attempt=round(mean(n_attempt),1),
                                 pct_1000=round(100*mean(n_attempt==1000),1)), by=dgp])

cat("\n########## E. why replicates are lost, PER METHOD (H0 arm, own set)\n")
D <- unique(R[set=="own", .(dgp,model,N,ntime,label,method,n_present,n_conv,n_usable,n_no_se,n_big,n_attempt)])
E <- D[label=="manual", .(cells=.N,
       attempts = sum(n_attempt),
       present = sum(n_present),
       nonconv = sum(n_present - n_conv),
       no_se   = sum(n_no_se),
       big_est = sum(n_big),
       usable  = sum(n_usable)), by=method]
E[, `:=`(pct_absent = round(100*(attempts-present)/attempts,2),
         pct_nonconv = round(100*nonconv/attempts,2),
         pct_no_se = round(100*no_se/attempts,2),
         pct_big = round(100*big_est/attempts,2),
         pct_usable = round(100*usable/attempts,2))]
print(as.data.frame(E[order(-pct_usable)]), row.names=FALSE)
cat("\nsame at N<=20 only:\n")
E2 <- D[label=="manual" & N<=20, .(attempts=sum(n_attempt), no_se=sum(n_no_se), big=sum(n_big),
        nonconv=sum(n_present-n_conv), absent=sum(n_attempt-n_present), usable=sum(n_usable)), by=method]
E2[, `:=`(pct_no_se=round(100*no_se/attempts,2), pct_big=round(100*big/attempts,2),
          pct_nonconv=round(100*nonconv/attempts,2), pct_absent=round(100*absent/attempts,2),
          pct_usable=round(100*usable/attempts,2))]
print(as.data.frame(E2[order(-pct_usable)]), row.names=FALSE)

cat("\n-- how many replicates does the |estimate|>=10 / se>=10 filter remove, per method? --\n")
cat("(n_big counts |estimate|>=10 among CONVERGED rows; n_no_se counts missing SE among converged)\n")
tot <- D[, .(attempts=sum(n_attempt), big=sum(n_big), no_se=sum(n_no_se)), by=method]
tot[, `:=`(pct_big=round(100*big/attempts,3), pct_no_se=round(100*no_se/attempts,3))]
print(as.data.frame(tot[order(-big)]), row.names=FALSE)
cat("\n-- the D-07 claim: 44,169 LRT rows a Wald SE filter would have discarded --\n")
lr <- R[set=="own" & method=="glmer_lrt", .(usable=sum(n_usable)), by=.(dgp)]
print(as.data.frame(lr), row.names=FALSE)
lrn <- unique(R[set=="own" & method=="glmer_lrt", .(dgp,model,N,ntime,label,n_no_se,n_big,n_usable)])
cat(sprintf("LRT rows with n_no_se>0: %d cells; total no_se flagged for LRT: %d (by construction 0)\n",
  sum(lrn$n_no_se>0), sum(lrn$n_no_se)))

cat("\n########## C. coverage of GLMM-native with the T-varying marginal truth\n")
gh <- gaussHermiteData(80)
marg_p <- function(eta, V) sum(gh$w*plogis(eta + sqrt(2*V)*gh$x))/sqrt(pi)
implied <- function(bint, T, b0, bt, btr, Dm, re) {
  g <- expand.grid(time=0:(T-1), trt=0:1); X <- model.matrix(~time*trt, g)
  Z <- if (re=="ri") matrix(1,nrow(g),1) else cbind(1,g$time)
  Dm <- as.matrix(Dm); V <- rowSums((Z%*%Dm)*Z)
  eta <- b0 + bt*g$time + btr*g$trt + bint*g$time*g$trt
  p <- vapply(seq_along(eta), function(i) marg_p(eta[i], V[i]), numeric(1))
  beta <- rep(0, ncol(X))
  for (it in 1:400) { mu <- plogis(X%*%beta); W <- as.vector(mu*(1-mu))
    z <- X%*%beta + (p-mu)/pmax(W,1e-12)
    nb <- solve(crossprod(X,X*W), crossprod(X,W*z)); if (max(abs(nb-beta))<1e-13){beta<-nb;break}; beta<-nb }
  as.vector(beta)[4]
}
D3 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
TL <- sort(unique(R$ntime))
TRUTH <- rbindlist(lapply(TL, function(T) rbindlist(lapply(1:5, function(m) {
  if (m %in% c(1,2)) tm <- implied(if(m==4) 0.058 else 0.232, T, -2.32,
                                  if (m==4) 0.232 else 0.116, 0.232, 1.0, "ri")
  else tm <- implied(if(m==4) 0.058 else 0.232, T, -2.32,
                     if (m==4) 0.232 else 0.116, 0.232, D3, "rs")
  data.table(model=m, ntime=T, truth_T=tm)
}))))
cat("marginal truth by scenario and T (GLMM-native, exact):\n")
print(as.data.frame(dcast(TRUTH, model~ntime, value.var="truth_T")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,5) else x)]), row.names=FALSE)
cat(sprintf("range over T for M1/M3: %.5f - %.5f ; flat value used by agg_all5.R: 0.200\n",
  min(TRUTH[model %in% c(1,3)]$truth_T), max(TRUTH[model %in% c(1,3)]$truth_T)))

# coverage cannot be recomputed from k_cov (that used a flat truth), so recompute
# from the raw chunks -- done in REV2M_e2.R. Here: bound the change analytically
# using mean_se, mean_est and sd_est per cell (normal approximation), plus the
# EXACT recomputation is in e2.
cat("\n(exact recomputation from the chunks: see REV2M_e2.R)\n")

cat("\n########## F. size-adjusted power: is the producing code on the server?\n")
fs <- system("grep -rln 'power_adj\\|size.adjust\\|sizeadj' ~/MMM --include='*.R' --include='*.Rmd' --include='*.sh' --include='*.md' 2>/dev/null | grep -v renv", intern=TRUE)
cat("files mentioning size-adjusted power:\n"); print(fs)
cat("\nP_power_sizeadj.csv exists:", file.exists(path.expand("~/MMM/tables/P_power_sizeadj.csv")), "\n")
pp <- fread(path.expand("~/MMM/tables/P_power_sizeadj.csv"))
cat("columns:", paste(names(pp), collapse=", "), "\n")
cat("its granularity: rows per dgp x method x N; no `model`/`ntime` column =>\n")
print(pp[, .(rows=.N, N_levels=uniqueN(N), methods=uniqueN(method)), by=dgp])
cat("\nrange of power_adj by dgp, pooled over N (replicate-weighted):\n")
print(as.data.frame(pp[, .(adj=round(sum(power_adj*n)/sum(n),2), raw=round(sum(power_raw*n)/sum(n),2)),
                       by=.(dgp,method)][order(dgp,-adj)]), row.names=FALSE)
cat("\nDONE_E\n")
