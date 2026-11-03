# REV2M_a4.R -- seeds, truncation, estimability recomputation (no MMLB needed).
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(data.table); library(MMLB); library(MASS)})
SRC <- path.expand("~/MMM/R")
for (m in c(1,3)) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=globalenv())
sys.source(file.path(SRC,"generate_data.R"), envir=globalenv())
sys.source(file.path(SRC,"make_all_data.R"), envir=globalenv())

cat("########## SEEDS: set.seed(20260805 + 1000*model + i)\n")
g <- function(N,T,m,lab,rep) { set.seed(20260805L+1000L*m+rep); make_all_data(N,T,rep,as.character(m),lab) }
d10 <- g(10,4,3,"mmlb",1); d15 <- g(15,4,3,"mmlb",1)
cat(sprintf("N=10 vs N=15 (m3,T4,rep1): trt of subjects 1-10 identical? %s | cluster sizes identical? %s\n",
  identical(unname(tapply(d10$trt,d10$id,`[`,1))[1:10], unname(tapply(d15$trt,d15$id,`[`,1))[1:10]),
  identical(unname(c(table(d10$id)))[1:10], unname(c(table(d15$id)))[1:10])))
t4 <- g(45,4,3,"mmlb",2); t10 <- g(45,10,3,"mmlb",2)
cat(sprintf("T=4 vs T=10 (same m,rep): trt vector over the 45 subjects identical? %s\n",
  identical(unname(tapply(t4$trt,t4$id,`[`,1)), unname(tapply(t10$trt,t10$id,`[`,1)))))
h1 <- g(45,4,3,"mmlb",7); h0 <- g(45,4,3,"manual",7)
cat(sprintf("H1 vs H0 (same cell, same rep): design (id,time,trt) identical? %s | outcome identical? %s\n",
  identical(h1[,c("id","time","trt")], h0[,c("id","time","trt")]), identical(h1$binY, h0$binY)))
cat(sprintf("H1/H0 prevalence in that single replicate: %.4f vs %.4f\n", mean(h1$binY), mean(h0$binY)))
cat(sprintf("cells sharing each seed value: 16 N x 6 T x 2 labels x 3 mechanisms = %d\n", 16*6*2*3))
cat(sprintf("seed ranges by model: m=1 [%d,%d], m=2 [%d,%d] -> overlap? %s\n",
  20260805+1000+1, 20260805+1000+1000, 20260805+2000+1, 20260805+2000+1000,
  (20260805+1000+1000) >= (20260805+2000+1)))

cat("\n########## TRUNCATION\n")
tr <- do.call(rbind, lapply(c(2,3,4,6,8,10), function(T) {
  ns <- c(); for (i in 1:300) { set.seed(20260805L+1000L*3+i)
    nc <- ifelse(runif(500)<0.7, T, sample(seq_len(T-1), 500, replace=TRUE)); ns <- c(ns,nc) }
  data.frame(ntime=T, mean_obs=mean(ns), pct_full=100*mean(ns==T), pct_singleton=100*mean(ns==1)) }))
print(as.data.frame(tr %>% mutate(across(everything(), ~round(.x,3)))), row.names=FALSE)

cat("\n########## ESTIMABILITY (§3.8, §3.9) from REV_cells.csv\n")
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv")); R[, model := as.integer(model)]
D <- unique(R[set=="own", .(dgp,model,N,ntime,label,method,n_attempt,n_usable)])
D[, re := ifelse(model %in% c(1,2), "RI (M1-M2)", "RS (M3-M5)")]
tabE <- D[, .(pct=round(100*sum(n_usable)/sum(n_attempt),2)), by=.(re,ntime,method)]
print(as.data.frame(dcast(tabE, re+ntime~method, value.var="pct")), row.names=FALSE)
cat("\n-- RS, T=2: raw counts --\n")
print(as.data.frame(D[re=="RS (M3-M5)" & ntime==2,
  .(attempts=sum(n_attempt), usable=sum(n_usable), pct=round(100*sum(n_usable)/sum(n_attempt),4)),
  by=method][order(-pct)]), row.names=FALSE)
cat("\n-- RS, T=2, per hypothesis arm --\n")
print(as.data.frame(D[re=="RS (M3-M5)" & ntime==2 & method %in% c("glmer","shrink","glmer_lrt"),
  .(attempts=sum(n_attempt), usable=sum(n_usable)), by=.(method,label)]), row.names=FALSE)
cat("\n-- pooled over all 5 scenarios at T=2 (the 27.4% claim) --\n")
print(as.data.frame(D[ntime==2 & method %in% c("glmer","glmer_lrt","shrink"),
  .(attempts=sum(n_attempt), usable=sum(n_usable), pct=round(100*sum(n_usable)/sum(n_attempt),2)),
  by=method]), row.names=FALSE)
cat("\n-- §3.9 at N<=20: usable %, and design points a method never fits --\n")
d20 <- D[N<=20]
print(as.data.frame(d20[, .(usable_pct=round(100*sum(n_usable)/sum(n_attempt),2),
  dp_total=uniqueN(paste(dgp,model,N,ntime)),
  dp_with_fit=uniqueN(paste(dgp,model,N,ntime)[n_usable>0])), by=method][order(-usable_pct)]), row.names=FALSE)
cat("\n-- LRT level at N<=20 (the 8.9-11.1% claim) --\n")
L <- R[set=="own" & label=="manual" & N<=20, .(n=sum(n), k=sum(k_rej)), by=.(method,dgp)]
L[, rate := round(100*k/n,2)]
print(as.data.frame(dcast(L, method~dgp, value.var="rate")), row.names=FALSE)
cat("\nDONE_A4\n")
