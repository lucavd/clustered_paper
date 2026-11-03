# R3_mmm.R -- GIRO 3 punto 1: MMM con ENTRAMBE le varianze (mod.cov, rob.cov).
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(data.table); library(MMLB); library(MASS); library(parallel); library(fastGHQuad)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
for (m in c(1,3)) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=globalenv())
sys.source(file.path(SRC,"generate_data.R"), envir=globalenv())
sys.source(file.path(SRC,"make_all_data.R"), envir=globalenv())
sys.source(file.path(RD,"gen_glmm_data.R"), envir=globalenv())
eM <- new.env(parent=globalenv())
for (f in c("manual/generation_model1.R","manual/generation_model3.R","manual/generate_data.R",
            "manual/make_all_data.R")) sys.source(file.path(RD,f), envir=eM)
mk_native <- get("make_all_data", envir=eM)

zq <- qnorm(0.975)
MARG <- c("1"=0.200,"3"=0.200)

fit1 <- function(dgp, m, N, T, i, label) {
  d <- tryCatch({
    if (dgp=="GenBinaryY") { set.seed(20260805L+1000L*m+i); as.data.frame(make_all_data(N,T,i,as.character(m),label)) }
    else { set.seed(20260805L+1000L*m+i); as.data.frame(mk_native(N,T,i,as.character(m),label)) }
  }, error=function(e) NULL)
  if (is.null(d)) return(NULL)
  lvf <- if (m==1) ~1 else ~time
  mod <- tryCatch(suppressWarnings(mm(binY ~ time*trt, lv.formula=lvf, t.formula=~1,
                                      data=d, id=id, q=20)), error=function(e) NULL)
  if (is.null(mod)) return(NULL)
  tb <- tryCatch(summary(mod)$mean.table, error=function(e) NULL)
  if (is.null(tb) || !"time:trt" %in% rownames(tb)) return(NULL)
  k <- which(rownames(tb)=="time:trt")
  est <- tb[k,1]
  se_mod <- tryCatch(sqrt(diag(mod$mod.cov))[k], error=function(e) NA_real_)
  se_rob <- tryCatch(sqrt(diag(mod$rob.cov))[k], error=function(e) NA_real_)
  data.table(dgp=dgp, model=m, N=N, ntime=T, label=label, rep=i, est=est,
             se_tab=tb[k,2], se_mod=as.numeric(se_mod), se_rob=as.numeric(se_rob))
}
grid <- CJ(dgp=c("GenBinaryY","GLMM-native"), model=c(1,3), N=c(10,25,45,105), ntime=c(4,10),
           label=c("manual","mmlb"))
NREP <- 300L
t0 <- Sys.time(); out <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k]
  if (as.numeric(difftime(Sys.time(),t0,units="mins")) > 170) { cat("CAP a k=",k,"\n"); break }
  rr <- mclapply(1:NREP, function(i) fit1(g$dgp,g$model,g$N,g$ntime,i,g$label), mc.cores=20)
  out[[k]] <- rbindlist(Filter(Negate(is.null), rr))
  cat(sprintf("[%s] %s m%d N=%d T=%d %s: %d fit (%.1f min)\n", format(Sys.time(),"%H:%M:%S"),
    g$dgp,g$model,g$N,g$ntime,g$label, nrow(out[[k]]),
    as.numeric(difftime(Sys.time(),t0,units="mins")))); flush.console()
}
D <- rbindlist(out)
saveRDS(D, path.expand("~/MMM/rev2/R3_mmm_raw.rds"))
cat("\n=== conferma: se_tab e' mod.cov, non rob.cov ===\n")
cat(sprintf("fit totali %d ; max |se_tab - se_mod| %.3e ; max |se_tab - se_rob| %.3e\n",
  nrow(D), max(abs(D$se_tab-D$se_mod),na.rm=TRUE), max(abs(D$se_tab-D$se_rob),na.rm=TRUE)))
cat(sprintf("rob.cov disponibile in %d fit su %d (%.1f%%)\n", sum(is.finite(D$se_rob)), nrow(D),
  100*mean(is.finite(D$se_rob))))
U <- D[is.finite(est) & is.finite(se_mod) & is.finite(se_rob) & abs(est)<10 &
       se_mod>0 & se_mod<10 & se_rob>0 & se_rob<10]
cat(sprintf("fit utilizzabili con ENTRAMBE le varianze: %d\n", nrow(U)))
U[, `:=`(z_mod=est/se_mod, z_rob=est/se_rob,
         rej_mod=abs(est/se_mod)>zq, rej_rob=abs(est/se_rob)>zq,
         t_marg=MARG[as.character(model)])]
U[, `:=`(cov_mod=abs(est-t_marg)<=zq*se_mod, cov_rob=abs(est-t_marg)<=zq*se_rob)]

cat("\n===== LIVELLO (braccio H0) con SE model-based vs robusta =====\n")
L <- U[label=="manual", .(n=.N, lev_mod=100*mean(rej_mod), lev_rob=100*mean(rej_rob)), by=.(dgp,model,N,ntime)]
L[, `:=`(shift=round(lev_rob-lev_mod,2), mcse=round(100*sqrt(0.05*0.95/n),2))]
print(as.data.frame(L[order(dgp,model,ntime,N), .(dgp,model,N,ntime,n,
  lev_mod=round(lev_mod,2), lev_rob=round(lev_rob,2), shift, mcse)]), row.names=FALSE)
cat("\n-- pooled --\n")
P <- U[label=="manual", .(n=.N, lev_mod=round(100*mean(rej_mod),2), lev_rob=round(100*mean(rej_rob),2)), by=dgp]
P[, shift := round(lev_rob-lev_mod,2)]
print(as.data.frame(P), row.names=FALSE)
cat(sprintf("pooled complessivo: model-based %.2f%% | robusta %.2f%% | shift %+.2f pp (n=%d)\n",
  100*mean(U[label=="manual"]$rej_mod), 100*mean(U[label=="manual"]$rej_rob),
  100*(mean(U[label=="manual"]$rej_rob)-mean(U[label=="manual"]$rej_mod)), nrow(U[label=="manual"])))

cat("\n===== WORST CASE: massimo scostamento dal 5% nominale, per cella =====\n")
W <- L[, .(worst_mod=max(abs(lev_mod-5)), worst_rob=max(abs(lev_rob-5)),
           cella_mod=paste0("N=",N[which.max(abs(lev_mod-5))],",T=",ntime[which.max(abs(lev_mod-5))]),
           cella_rob=paste0("N=",N[which.max(abs(lev_rob-5))],",T=",ntime[which.max(abs(lev_rob-5))])), by=dgp]
print(as.data.frame(W[, .(dgp, worst_mod=round(worst_mod,2), cella_mod, worst_rob=round(worst_rob,2), cella_rob)]),
      row.names=FALSE)
cat(sprintf("worst case complessivo: model-based %.2f pp | robusta %.2f pp\n",
  max(abs(L$lev_mod-5)), max(abs(L$lev_rob-5))))
cat("-- worst case ai pochi cluster (N<=25), dove §3.1 assegna il primato --\n")
L25 <- L[N<=25]
cat(sprintf("N<=25: model-based %.2f pp | robusta %.2f pp\n", max(abs(L25$lev_mod-5)), max(abs(L25$lev_rob-5))))

cat("\n===== COPERTURA con le due varianze =====\n")
C <- U[label=="mmlb", .(n=.N, cov_mod=round(100*mean(cov_mod),2), cov_rob=round(100*mean(cov_rob),2)),
       by=.(dgp,model,N,ntime)]
C[, shift := round(cov_rob-cov_mod,2)]
print(as.data.frame(C[order(dgp,model,ntime,N)]), row.names=FALSE)
cat("\n-- pooled --\n")
print(as.data.frame(U[label=="mmlb", .(n=.N, cov_mod=round(100*mean(cov_mod),2),
  cov_rob=round(100*mean(cov_rob),2), shift=round(100*(mean(cov_rob)-mean(cov_mod)),2)), by=dgp]), row.names=FALSE)

cat("\n===== SE / SD empirica con le due varianze, per cella =====\n")
R <- U[, .(n=.N, sd_emp=sd(est), se_mod=mean(se_mod), se_rob=mean(se_rob)), by=.(dgp,model,N,ntime,label)]
R[, `:=`(r_mod=round(se_mod/sd_emp,4), r_rob=round(se_rob/sd_emp,4))]
print(as.data.frame(R[label=="manual"][order(dgp,model,ntime,N),
  .(dgp,model,N,ntime,n,sd_emp=round(sd_emp,4),se_mod=round(se_mod,4),se_rob=round(se_rob,4),r_mod,r_rob)]),
  row.names=FALSE)
cat(sprintf("\nrapporto SE/SD mediano: model-based %.4f | robusta %.4f\n",
  median(R$r_mod), median(R$r_rob)))
cat(sprintf("ai pochi cluster (N<=25): model-based %.4f | robusta %.4f\n",
  median(R[N<=25]$r_mod), median(R[N<=25]$r_rob)))
cat(sprintf("a molti cluster (N>=45): model-based %.4f | robusta %.4f\n",
  median(R[N>=45]$r_mod), median(R[N>=45]$r_rob)))
fwrite(L, path.expand("~/MMM/rev2/R3_mmm_level.csv"))
fwrite(C, path.expand("~/MMM/rev2/R3_mmm_cov.csv"))
fwrite(R, path.expand("~/MMM/rev2/R3_mmm_seratio.csv"))
cat("\nDONE_MMM3\n")
