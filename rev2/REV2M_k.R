# REV2M_k.R -- seed reuse across mechanisms: are the covariate designs identical?
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(MMLB); library(MASS); library(fastGHQuad)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
eG <- new.env(); eM <- new.env(); eA <- new.env()
for (m in 1:5) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=eG)
sys.source(file.path(SRC,"generate_data.R"), envir=eG); sys.source(file.path(SRC,"make_all_data.R"), envir=eG)
sys.source(file.path(RD,"gen_glmm_data.R"), envir=eM)
for (f in c("manual/generation_model1.R","manual/generation_model3.R","manual/generation_model5.R",
            "manual/generate_data.R","manual/make_all_data.R")) sys.source(file.path(RD,f), envir=eM)
sys.source(file.path(RD,"gen_ar1_data.R"), envir=eA)
for (f in c("ar1/generation_models_ar1.R","manual/generate_data.R","ar1/make_all_data.R"))
  sys.source(file.path(RD,f), envir=eA)
gG <- get("make_all_data", envir=eG); gM <- get("make_all_data", envir=eM); gA <- get("make_all_data", envir=eA)
for (m in c(1,3)) for (T in c(4,10)) {
  set.seed(20260805L+1000L*m+5); a <- gG(45,T,5,as.character(m),"mmlb")
  set.seed(20260805L+1000L*m+5); b <- gM(45,T,5,as.character(m),"mmlb")
  set.seed(20260805L+1000L*m+5); d <- gA(45,T,5,as.character(m),"mmlb")
  cat(sprintf("m=%d T=%2d | nrow G/M/A: %d/%d/%d | (id,time,trt) G==M: %s  G==A: %s  M==A: %s | outcome G==M: %s\n",
    m, T, nrow(a), nrow(b), nrow(d),
    identical(a[,c("id","time","trt")], b[,c("id","time","trt")]),
    identical(a[,c("id","time","trt")], as.data.frame(d)[,c("id","time","trt")]),
    identical(b[,c("id","time","trt")], as.data.frame(d)[,c("id","time","trt")]),
    identical(a$binY, b$binY)))
}
cat("\n-- are the three mechanisms' outcomes correlated within a replicate? --\n")
for (m in c(1,3)) {
  aa <- bb <- dd <- c()
  for (rep in 1:60) {
    set.seed(20260805L+1000L*m+rep); a <- gG(45,6,rep,as.character(m),"manual")
    set.seed(20260805L+1000L*m+rep); b <- gM(45,6,rep,as.character(m),"manual")
    set.seed(20260805L+1000L*m+rep); d <- gA(45,6,rep,as.character(m),"manual")
    n <- min(nrow(a), nrow(b), nrow(d))
    aa <- c(aa, mean(a$binY)); bb <- c(bb, mean(b$binY)); dd <- c(dd, mean(d$binY))
  }
  cat(sprintf("m=%d prevalence across 60 replicates: cor(G,M)=%.3f cor(G,A)=%.3f cor(M,A)=%.3f\n",
    m, cor(aa,bb), cor(aa,dd), cor(bb,dd)))
}
cat("\n-- design reuse: same (model, rep) across the 16 N and 6 T levels --\n")
set.seed(20260805L+1000L*3+9); x1 <- gG(10,4,9,"3","mmlb")
set.seed(20260805L+1000L*3+9); x2 <- gG(10,6,9,"3","mmlb")
cat(sprintf("N=10, T=4 vs T=6: first-subject trt %d vs %d ; cluster sizes %s vs %s\n",
  x1$trt[1], x2$trt[1], paste(head(unname(c(table(x1$id))),5), collapse=","),
  paste(head(unname(c(table(x2$id))),5), collapse=",")))
cat("\nDONE_K\n")
