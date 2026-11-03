# REV2M_a3.R -- A quantified by PATCHING the real generator: run the mmlb branch
# (H1 dependence) with the manual branch's beta vector. Everything else identical.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr); library(data.table)
                  library(MMLB); library(MASS); library(geepack); library(parallel)})
SRC <- path.expand("~/MMM/R")
for (m in 1:5) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=globalenv())
sys.source(file.path(SRC,"generate_data.R"), envir=globalenv())
sys.source(file.path(SRC,"make_all_data.R"), envir=globalenv())

# build a "hybrid" generator: H0 mean vector, H1 sigma/gamma
mk_hybrid <- function(m) {
  src <- deparse(get(sprintf("generate_data_model%d", m)))
  h0 <- switch(as.character(m),
    "3" = "c(-2, 0.1, 0.2, 0)", "4" = "c(-2, 0.2, 0.2, 0)", "5" = "c(-2, 0.1, 0.2, 0.05, 0)")
  h1 <- switch(as.character(m),
    "3" = "c(-2, 0.1, 0.2, 0.2)", "4" = "c(-2, 0.2, 0.2, 0.05)", "5" = "c(-2, 0.1, 0.2, 0.05, 0.2)")
  src2 <- sub(paste0("mmlb = ", h1), paste0("mmlb = ", h0), paste(src, collapse="\n"), fixed=TRUE)
  stopifnot(grepl(paste0("mmlb = ", h0), src2, fixed=TRUE))
  eval(parse(text=src2))
}
HY <- lapply(3:5, mk_hybrid); names(HY) <- as.character(3:5)
cat("hybrid generators built for M3-M5 (H0 beta, H1 sigma/gamma):",
    paste(names(HY), collapse=","), "\n")

dep <- function(d) {
  d <- d[order(d$id, d$time), ]
  s <- split(seq_len(nrow(d)), d$id); a<-b<-cc<-dd<-0
  for (ix in s) { y <- d$binY[ix]; if (length(y)<2) next
    for (j in 2:length(y)) { if (y[j-1]==1&&y[j]==1) a<-a+1 else if (y[j-1]==1&&y[j]==0) b<-b+1
      else if (y[j-1]==0&&y[j]==1) cc<-cc+1 else dd<-dd+1 } }
  c(prev=mean(d$binY), lag1_or=(a*dd)/max(b*cc,1))
}
gen <- function(m, T, mode, nrep=20, N=600) {
  acc <- vector("list", nrep)
  for (rep in 1:nrep) {
    set.seed(20260805L + 1000L*m + rep)
    acc[[rep]] <- if (mode=="hybrid") HY[[as.character(m)]](N, T, rep, "mmlb")
                  else make_all_data(N, T, rep, as.character(m), mode)
  }
  bind_rows(acc)
}
rows <- list(); k<-0
for (m in 3:5) for (T in c(4,10)) {
  a1 <- dep(gen(m,T,"mmlb")); a0 <- dep(gen(m,T,"manual")); ah <- dep(gen(m,T,"hybrid"))
  k<-k+1; rows[[k]] <- data.frame(model=m, ntime=T,
    prev_H1=a1["prev"], prev_H0=a0["prev"], prev_H0dep1=ah["prev"],
    or_H1=a1["lag1_or"], or_H0=a0["lag1_or"], or_H0dep1=ah["lag1_or"])
}
A <- bind_rows(rows)
cat("\n===== prevalence and lag-1 odds ratio: H1 | H0 as run | H0 mean with H1 dependence =====\n")
print(as.data.frame(A %>% mutate(across(where(is.numeric), ~round(.x,4)))), row.names=FALSE)
cat("\n===== the part of the H1-vs-H0 gap that is NOT the tested parameter =====\n")
B <- A %>% transmute(model, ntime,
  gap_or_total = round(or_H1 - or_H0, 4),
  gap_or_from_dependence = round(or_H0dep1 - or_H0, 4),
  share_dependence_pct = round(100*(or_H0dep1-or_H0)/(or_H1-or_H0), 1),
  gap_prev_total = round(prev_H1 - prev_H0, 4),
  gap_prev_from_dependence = round(prev_H0dep1 - prev_H0, 4))
print(as.data.frame(B), row.names=FALSE)
write.csv(A, path.expand("~/MMM/rev2/REV2M_arm_dep.csv"), row.names=FALSE)

cat("\n===== consequence for the LEVEL: GEE type I error under H0-as-run vs H0-with-H1-dependence =====\n")
lev <- function(m, T, N, mode, nrep=400) {
  ps <- unlist(mclapply(1:nrep, function(rep) {
    set.seed(20260805L + 1000L*m + rep)
    d <- tryCatch(if (mode=="hybrid") HY[[as.character(m)]](N,T,rep,"mmlb")
                  else make_all_data(N,T,rep,as.character(m),mode), error=function(e) NULL)
    if (is.null(d)) return(NA_real_)
    d <- as.data.frame(d); d$id <- factor(d$id)
    f <- if (m==5) binY ~ time*trt + dose else binY ~ time*trt
    fit <- tryCatch(geeglm(f, id=id, data=d, family=binomial, corstr="exchangeable", scale.fix=TRUE),
                    error=function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    s <- summary(fit)$coefficients
    if (!"time:trt" %in% rownames(s)) return(NA_real_)
    est <- s["time:trt",1]; se <- s["time:trt",2]
    if (!is.finite(est)||!is.finite(se)||abs(est)>=10||se<=0||se>=10) return(NA_real_)
    s["time:trt",4]
  }, mc.cores=20))
  ps <- ps[is.finite(ps)]
  c(n=length(ps), rate=100*mean(ps<0.05))
}
LV <- list(); k<-0
for (m in c(3,5)) for (T in c(4,10)) for (N in c(20,45)) {
  a <- lev(m,T,N,"manual"); b <- lev(m,T,N,"hybrid")
  k<-k+1; LV[[k]] <- data.frame(model=m,ntime=T,N=N,
    level_asrun=a["rate"], n_asrun=a["n"], level_dep1=b["rate"], n_dep1=b["n"],
    shift_pp=round(b["rate"]-a["rate"],2))
}
L <- bind_rows(LV)
print(as.data.frame(L %>% mutate(across(where(is.numeric), ~round(.x,2)))), row.names=FALSE)
cat(sprintf("\nmean |shift| %.2f pp, max %.2f pp over %d probed cells\n",
            mean(abs(L$shift_pp)), max(abs(L$shift_pp)), nrow(L)))
write.csv(L, path.expand("~/MMM/rev2/REV2M_arm_level.csv"), row.names=FALSE)
cat("\nDONE_A3\n")
