# R2_mmlb.R -- controllo 7: quale covarianza finisce nella mean.table di mm()?
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(MMLB); library(dplyr); library(data.table)})
SRC <- path.expand("~/MMM/R")
for (f in c("generation_model3.R","generate_data.R","make_all_data.R"))
  sys.source(file.path(SRC,f), envir=globalenv())

cat("=== summary() per la classe di mm(): quale slot usa? ===\n")
set.seed(20260805L+3000L+1); d <- as.data.frame(make_all_data(60, 6, 1, "3", "mmlb"))
mod <- mm(binY ~ time*trt, lv.formula=~time, t.formula=~1, data=d, id=id, q=20)
cat("classe:", class(mod), "\n")
cat("componenti dell'oggetto:\n"); print(names(mod))
s <- summary(mod)
cat("componenti di summary():\n"); print(names(s))
mt <- s$mean.table
cat("\nmean.table:\n"); print(round(mt, 6))
se_pub <- mt["time:trt", 2]
cat("\nSE della mean.table per time:trt:", format(se_pub, digits=10), "\n")
for (nm in intersect(c("mod.cov","modelcov","cov.model","empiricalcov","cov.empirical",
                       "robust.cov","rob.cov","Cheese","cov"), names(mod))) {
  M <- mod[[nm]]
  if (is.matrix(M)) {
    k <- which(rownames(mt)=="time:trt")
    v <- tryCatch(sqrt(diag(M))[k], error=function(e) NA)
    cat(sprintf("  slot '%s': %dx%d, sqrt(diag)[%d] = %s  -> uguale alla mean.table? %s\n",
      nm, nrow(M), ncol(M), k, format(v, digits=10),
      isTRUE(all.equal(as.numeric(v), as.numeric(se_pub), tolerance=1e-8))))
  }
}
cat("\n=== confronto numerico model-based vs empirica (sandwich) sulla stessa stima ===\n")
for (nm in intersect(c("mod.cov","modelcov","empiricalcov","cov.empirical","Cheese"), names(mod))) {
  M <- mod[[nm]]
  if (is.matrix(M)) {
    k <- which(rownames(mt)=="time:trt")
    cat(sprintf("  %s: SE = %.6f\n", nm, sqrt(diag(M))[k]))
  }
}
cat("\n=== verifica empirica su 400 repliche: SE media di mm() vs SD empirica ===\n")
suppressMessages(library(parallel))
rr <- unlist(mclapply(1:400, function(i) {
  set.seed(20260805L+3000L+i)
  dd <- tryCatch(as.data.frame(make_all_data(60, 6, i, "3", "mmlb")), error=function(e) NULL)
  if (is.null(dd)) return(NA)
  m <- tryCatch(mm(binY ~ time*trt, lv.formula=~time, t.formula=~1, data=dd, id=id, q=20),
                error=function(e) NULL)
  if (is.null(m)) return(c(NA,NA))
  t <- summary(m)$mean.table
  if (!"time:trt" %in% rownames(t)) return(c(NA,NA))
  c(t["time:trt",1], t["time:trt",2])
}, mc.cores=20))
M <- matrix(rr, ncol=2, byrow=TRUE); M <- M[is.finite(M[,1]) & is.finite(M[,2]) & abs(M[,1])<10 & M[,2]<10, ]
cat(sprintf("repliche utilizzabili %d ; SE media %.5f ; SD empirica %.5f ; rapporto %.4f\n",
  nrow(M), mean(M[,2]), sd(M[,1]), mean(M[,2])/sd(M[,1])))
cat("(un rapporto vicino a 1 e' compatibile con una varianza correttamente calibrata su questo DGP,\n")
cat(" che qui non e' misspecificato: MMLB fitta il modello generatore.)\n")
cat("\nDONE_MMLB\n")
