# REV2M_a2.R -- A (completed), plus seeds, truncation, estimability table.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr); library(data.table)
                  library(MMLB); library(MASS); library(geepack); library(parallel)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
for (m in 1:5) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=globalenv())
sys.source(file.path(SRC,"generate_data.R"), envir=globalenv())
sys.source(file.path(SRC,"make_all_data.R"), envir=globalenv())

cat("########## A. GenBinaryY: sigma/gamma by arm, all five scenarios\n")
tabS <- data.frame(
  model = 1:5,
  sigma_H1 = c("1","1","c(1, 0.01)","c(1, 0.01)","c(1, 0.01, 0.01)"),
  sigma_H0 = c("1","1","c(1, 0)","c(1, 0)","c(1, 0, 0)"),
  gamma_H1 = c("1","1","1","c(1, 0.01)","c(1, 0.1, 0.05)"),
  gamma_H0 = c("1","1","1","c(1, 0)","c(1, 0, 0)"),
  identical_arms = c(TRUE,TRUE,FALSE,FALSE,FALSE))
print(tabS, row.names=FALSE)

dep <- function(d) {
  d <- d[order(d$id, d$time), ]
  s <- split(seq_len(nrow(d)), d$id); a<-b<-cc<-dd<-0
  for (ix in s) { y <- d$binY[ix]; if (length(y)<2) next
    for (j in 2:length(y)) { if (y[j-1]==1&&y[j]==1) a<-a+1 else if (y[j-1]==1&&y[j]==0) b<-b+1
      else if (y[j-1]==0&&y[j]==1) cc<-cc+1 else dd<-dd+1 } }
  c(prev=mean(d$binY), lag1_or=(a*dd)/max(b*cc,1))
}
gen_many <- function(m, T, lab, sigma_override=NULL, gamma_override=NULL, nrep=20, N=600) {
  acc <- vector("list", nrep)
  for (rep in 1:nrep) {
    set.seed(20260805L + 1000L*m + rep)
    if (is.null(sigma_override)) {
      acc[[rep]] <- make_all_data(N, T, rep, as.character(m), lab)
    } else {
      # H0 mean vector but H1 dependence: isolates the dependence asymmetry
      bv <- switch(as.character(m),
        "3"=c(-2,0.1,0.2,0), "4"=c(-2,0.2,0.2,0), "5"=c(-2,0.1,0.2,0.05,0))
      nclust <- ifelse(runif(N) < 0.7, T, sample(seq_len(T-1), N, replace=TRUE))
      id <- rep(seq(N), nclust); trt <- rep(rbinom(N,1,0.5), nclust)
      time <- unlist(lapply(nclust, function(x) seq(x)-1))
      age <- scale(rep(rnorm(N,50,10), nclust))
      dat <- data.frame(id, time, trt, age)
      if (m==5) { dat$sex <- rep(rbinom(N,1,0.5), nclust)
                  dat$steroids <- rep(sample(1:3,N,replace=TRUE), nclust)
                  dat$dose <- as.numeric(scale(rnorm(length(id)))) }
      dat <- dat %>% arrange(id, time)
      mf <- if (m==5) ~ time*trt + dose else ~ time*trt
      lv <- if (m==5) ~ time + dose else ~ time
      tf <- switch(as.character(m), "3"=~1, "4"=~trt, "5"=~trt+steroids)
      acc[[rep]] <- GenBinaryY(mean.formula=mf, lv.formula=lv, t.formula=tf,
        beta=bv, sigma=sigma_override, gamma=gamma_override,
        id=dat$id, data=dat, q=20, Yname="binY")
    }
  }
  bind_rows(acc)
}
rows <- list(); k <- 0
for (m in 3:5) for (T in c(4,10)) {
  s1 <- switch(as.character(m), "3"=c(1,0.01), "4"=c(1,0.01), "5"=c(1,0.01,0.01))
  g1 <- switch(as.character(m), "3"=1,         "4"=c(1,0.01), "5"=c(1,0.1,0.05))
  d1 <- dep(gen_many(m,T,"mmlb"))
  d0 <- dep(gen_many(m,T,"manual"))
  dM <- dep(gen_many(m,T,"manual", sigma_override=s1, gamma_override=g1))
  k<-k+1; rows[[k]] <- data.frame(model=m, ntime=T,
    prev_H1=d1["prev"], prev_H0_asrun=d0["prev"], prev_H0_matched=dM["prev"],
    or_H1=d1["lag1_or"], or_H0_asrun=d0["lag1_or"], or_H0_matched=dM["lag1_or"])
}
A <- bind_rows(rows)
cat("\n-- prevalence and lag-1 odds ratio: H1, H0 as run, H0 with H1's dependence --\n")
print(as.data.frame(A %>% mutate(across(where(is.numeric), ~round(.x,4)))), row.names=FALSE)
cat("\n-- the asymmetry attributable to sigma/gamma alone (H0 matched minus H0 as run) --\n")
print(as.data.frame(A %>% transmute(model, ntime,
  d_or_from_dependence = round(or_H0_matched - or_H0_asrun, 4),
  d_or_H1_vs_H0asrun   = round(or_H1 - or_H0_asrun, 4),
  share_pct = round(100*(or_H0_matched-or_H0_asrun)/(or_H1-or_H0_asrun),1))), row.names=FALSE)

cat("\n-- consequence for the level: GEE type I error, H0 as run vs H0 with H1's dependence --\n")
lev_one <- function(m, T, N, matched, nrep=300) {
  s1 <- switch(as.character(m), "3"=c(1,0.01), "4"=c(1,0.01), "5"=c(1,0.01,0.01))
  g1 <- switch(as.character(m), "3"=1,         "4"=c(1,0.01), "5"=c(1,0.1,0.05))
  ps <- mclapply(1:nrep, function(rep) {
    set.seed(20260805L + 1000L*m + rep)
    d <- tryCatch({
      if (!matched) make_all_data(N,T,rep,as.character(m),"manual") else {
        bv <- switch(as.character(m), "3"=c(-2,0.1,0.2,0), "4"=c(-2,0.2,0.2,0), "5"=c(-2,0.1,0.2,0.05,0))
        nclust <- ifelse(runif(N)<0.7, T, sample(seq_len(T-1), N, replace=TRUE))
        id <- rep(seq(N),nclust); trt <- rep(rbinom(N,1,0.5),nclust)
        time <- unlist(lapply(nclust, function(x) seq(x)-1))
        age <- scale(rep(rnorm(N,50,10),nclust)); dat <- data.frame(id,time,trt,age)
        if (m==5) { dat$sex <- rep(rbinom(N,1,0.5),nclust)
                    dat$steroids <- rep(sample(1:3,N,replace=TRUE),nclust)
                    dat$dose <- as.numeric(scale(rnorm(length(id)))) }
        dat <- dat %>% arrange(id,time)
        mf <- if (m==5) ~ time*trt + dose else ~ time*trt
        lv <- if (m==5) ~ time + dose else ~ time
        tf <- switch(as.character(m),"3"=~1,"4"=~trt,"5"=~trt+steroids)
        GenBinaryY(mean.formula=mf, lv.formula=lv, t.formula=tf, beta=bv,
                   sigma=s1, gamma=g1, id=dat$id, data=dat, q=20, Yname="binY") }
    }, error=function(e) NULL)
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
  }, mc.cores=20)
  p <- unlist(ps); p <- p[is.finite(p)]
  c(n=length(p), rate=100*mean(p<0.05))
}
LV <- list(); k<-0
for (m in c(3,5)) for (T in c(4,10)) for (N in c(20,45)) {
  a <- lev_one(m,T,N,FALSE); b <- lev_one(m,T,N,TRUE)
  k<-k+1; LV[[k]] <- data.frame(model=m, ntime=T, N=N,
     level_H0_asrun=a["rate"], n_asrun=a["n"],
     level_H0_matched=b["rate"], n_matched=b["n"])
}
L <- bind_rows(LV)
L$diff_pp <- round(L$level_H0_matched - L$level_H0_asrun, 2)
print(as.data.frame(L %>% mutate(across(where(is.numeric), ~round(.x,2)))), row.names=FALSE)
cat(sprintf("\nmean |shift| in GEE level from the arm asymmetry: %.2f pp (max %.2f)\n",
            mean(abs(L$diff_pp)), max(abs(L$diff_pp))))
write.csv(L, path.expand("~/MMM/rev2/REV2M_arm_level.csv"), row.names=FALSE)

cat("\n########## SEEDS\n")
d10 <- {set.seed(20260805L+1000L*3+1); make_all_data(10,4,1,"3","mmlb")}
d15 <- {set.seed(20260805L+1000L*3+1); make_all_data(15,4,1,"3","mmlb")}
cat(sprintf("N=10 vs N=15 (m3,T4,rep1): trt of subjects 1-10 identical? %s ; cluster sizes identical? %s\n",
  identical(unname(tapply(d10$trt,d10$id,`[`,1))[1:10], unname(tapply(d15$trt,d15$id,`[`,1))[1:10]),
  identical(unname(table(d10$id))[1:10], unname(table(d15$id))[1:10])))
dT4 <- {set.seed(20260805L+1000L*3+2); make_all_data(45,4,2,"3","mmlb")}
dT10 <- {set.seed(20260805L+1000L*3+2); make_all_data(45,10,2,"3","mmlb")}
cat(sprintf("T=4 vs T=10 same rep: trt vector of subjects 1-45 identical? %s\n",
  identical(unname(tapply(dT4$trt,dT4$id,`[`,1)), unname(tapply(dT10$trt,dT10$id,`[`,1)))))
h1 <- {set.seed(20260805L+1000L*3+7); make_all_data(45,4,7,"3","mmlb")}
h0 <- {set.seed(20260805L+1000L*3+7); make_all_data(45,4,7,"3","manual")}
cat(sprintf("H1 vs H0 same cell: design identical? %s ; outcome identical? %s\n",
  identical(h1[,c("id","time","trt")], h0[,c("id","time","trt")]), identical(h1$binY, h0$binY)))
cat(sprintf("cells sharing one seed value (same model, same replicate index): 16 N x 6 T x 2 labels x 3 dgp = %d\n", 16*6*2*3))
cat(sprintf("seed collision ACROSS models? model m spans 20260805+1000m+[1,1000]; m and m+1 overlap? %s\n",
  (20260805+1000*1+1000) >= (20260805+1000*2+1)))

cat("\n########## TRUNCATION (30%% of clusters lose follow-up)\n")
tr <- do.call(rbind, lapply(c(2,3,4,6,8,10), function(T) {
  ns <- c(); for (i in 1:300) { set.seed(20260805L+1000L*3+i)
    nc <- ifelse(runif(500)<0.7, T, sample(seq_len(T-1), 500, replace=TRUE)); ns <- c(ns,nc) }
  data.frame(ntime=T, mean_obs=mean(ns), pct_full=100*mean(ns==T),
             pct_singleton=100*mean(ns==1), pct_lt2=100*mean(ns<2)) }))
print(as.data.frame(tr %>% mutate(across(everything(), ~round(.x,3)))), row.names=FALSE)

cat("\n########## ESTIMABILITY table (§3.8/3.9) recomputed from REV_cells.csv\n")
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
R[, model := as.integer(model)]
D <- unique(R[set=="own", .(dgp,model,N,ntime,label,method,n_attempt,n_usable)])
D[, re := ifelse(model %in% c(1,2), "M1-M2 (RI)", "M3-M5 (RS)")]
tabE <- D[, .(usable_pct = round(100*sum(n_usable)/sum(n_attempt),2),
              attempts = sum(n_attempt), usable = sum(n_usable)), by=.(re, ntime, method)]
cat("-- usable %, pooled over mechanisms / cluster counts / hypotheses --\n")
print(as.data.frame(dcast(tabE, re+ntime~method, value.var="usable_pct")), row.names=FALSE)
cat("\n-- glmer and shrink at 2 occasions with random slopes: raw counts --\n")
print(as.data.frame(D[re=="M3-M5 (RS)" & ntime==2 & method %in% c("glmer","shrink","glmer_lrt","mmm","gee_fg"),
  .(attempts=sum(n_attempt), usable=sum(n_usable), pct=round(100*sum(n_usable)/sum(n_attempt),4)), by=method]),
  row.names=FALSE)
cat("\n-- same, H0 arm only (144,000 = 3 models x 3 dgp x 16 N x 1000) --\n")
print(as.data.frame(D[re=="M3-M5 (RS)" & ntime==2 & label=="manual" &
  method %in% c("glmer","shrink","glmer_lrt"),
  .(attempts=sum(n_attempt), usable=sum(n_usable)), by=method]), row.names=FALSE)
cat("\n-- pooled over all five scenarios at 2 occasions (the 27.4% claim) --\n")
print(as.data.frame(D[ntime==2 & method %in% c("glmer","glmer_lrt"),
  .(attempts=sum(n_attempt), usable=sum(n_usable), pct=round(100*sum(n_usable)/sum(n_attempt),2)),
  by=method]), row.names=FALSE)
cat("\n-- §3.9: N<=20, all 270 design points, cells never fitted counted at zero --\n")
d20 <- D[N<=20]
print(as.data.frame(d20[, .(usable_pct=round(100*sum(n_usable)/sum(n_attempt),2),
  design_points_with_any=uniqueN(paste(dgp,model,N,ntime)[n_usable>0]),
  design_points=uniqueN(paste(dgp,model,N,ntime))), by=method][order(-usable_pct)]), row.names=FALSE)
cat("\nDONE_A2\n")
