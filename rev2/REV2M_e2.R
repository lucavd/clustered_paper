# REV2M_e2.R -- C exact: recompute GLMM-native and AR1 coverage from the raw
# chunks under (i) the flat truth agg_all5.R used, (ii) the T-varying marginal
# truth for GLMM-native, (iii) the T-varying conditional truth for the GLMM arm.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(data.table); library(parallel); library(fastGHQuad)})
RES <- path.expand("~/MMM/results_full")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
zq <- qnorm(0.975)
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)
COND <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232)

# ---- T-varying truths -------------------------------------------------------
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
TL <- c(2,3,4,6,8,10)
TR <- rbindlist(lapply(TL, function(T) rbindlist(lapply(1:5, function(m) {
  bint <- if (m==4) 0.058 else 0.232
  bt   <- if (m==4) 0.232 else 0.116
  re   <- if (m %in% c(1,2)) "ri" else "rs"
  Dm   <- if (m %in% c(1,2)) 1.0 else D3
  data.table(model=m, ntime=T, truth_T=implied(bint, T, -2.32, bt, 0.232, Dm, re))
}))))
fwrite(TR, path.expand("~/MMM/rev2/REV2M_truth_byT.csv"))

do_dgp <- function(pref) {
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  one <- function(f) {
    x <- try(as.data.table(readRDS(f)), silent=TRUE)
    if (inherits(x,"try-error")) return(NULL)
    x <- x[term=="time:trt" & label=="mmlb" & method %in% WALD]
    if (!nrow(x)) return(NULL)
    x[, model := as.integer(model)]
    x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
    x <- x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10 & is.finite(p)]
    if (!nrow(x)) return(NULL)
    x[, .(method,model,N,ntime,estimate,se)]
  }
  L <- mclapply(fs, one, mc.cores=45)
  d <- rbindlist(L[!sapply(L,is.null)], fill=TRUE)
  d <- merge(d, TR, by=c("model","ntime"), all.x=TRUE)
  d[, t_flat_marg := MARG[as.character(model)]]
  d[, t_flat_cond := COND[as.character(model)]]
  d[, c_flat := ifelse(method=="glmer", abs(estimate-t_flat_cond)<=zq*se,
                                        abs(estimate-t_flat_marg)<=zq*se)]
  d[, c_flatmarg_all := abs(estimate-t_flat_marg)<=zq*se]
  d[, c_Tmarg := abs(estimate-truth_T)<=zq*se]   # T-varying marginal truth
  list(
    pooled = d[, .(n=.N,
      cov_flat = 100*mean(c_flat),
      cov_flat_marg = 100*mean(c_flatmarg_all),
      cov_Tmarg = 100*mean(c_Tmarg)), by=method][order(-cov_flat)],
    byT = d[, .(n=.N, cov_flat=100*mean(c_flat), cov_Tmarg=100*mean(c_Tmarg)), by=.(method,ntime)],
    byM = d[, .(n=.N, cov_flat=100*mean(c_flat), cov_Tmarg=100*mean(c_Tmarg)), by=.(method,model)],
    rmse = d[, .(n=.N,
      rmse_flat = sqrt(mean((estimate - ifelse(method=="glmer", t_flat_cond, t_flat_marg))^2)),
      rmse_Tmarg = sqrt(mean((estimate - truth_T)^2)),
      bias_flat = mean(estimate - ifelse(method=="glmer", t_flat_cond, t_flat_marg)),
      bias_Tmarg = mean(estimate - truth_T)), by=method]
  )
}

for (pref in c("GLMMnative","GLMMAR1","GenBinaryY")) {
  cat("\n##############################", pref, "##############################\n")
  o <- do_dgp(pref)
  cat("--- pooled coverage: flat truth (as published) vs T-varying marginal truth ---\n")
  print(as.data.frame(o$pooled[, lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)]), row.names=FALSE)
  cat("--- by T ---\n")
  print(as.data.frame(dcast(o$byT, method~ntime, value.var=c("cov_flat","cov_Tmarg"))[,
    lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)]), row.names=FALSE)
  cat("--- RMSE and bias, both targets ---\n")
  print(as.data.frame(o$rmse[, lapply(.SD, function(x) if(is.numeric(x)) round(x,4) else x)]), row.names=FALSE)
  saveRDS(o, path.expand(paste0("~/MMM/rev2/REV2M_cov_", pref, ".rds")))
}
cat("\nDONE_E2\n")
