# REV2M_a.R -- METHODS track, review nr.2. Cheap checks: no model fitting.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(MMLB); library(MASS); library(fastGHQuad)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
sp <- function(p) sys.source(p, envir = globalenv())

cat("################ A. GenBinaryY: does the H0 arm differ ONLY in the tested beta?\n")
for (m in 3:5) sp(file.path(SRC, sprintf("generation_model%d.R", m)))
sp(file.path(SRC, "generate_data.R")); sp(file.path(SRC, "make_all_data.R"))

# 1. declared parameter vectors, read straight out of the sourced functions
for (m in 3:5) {
  f <- get(sprintf("generate_data_model%d", m))
  b <- body(f)
  cat(sprintf("\n--- GenBinaryY M%d, switch() branches as sourced ---\n", m))
  txt <- deparse(b)
  ix <- grep("sigma <-|gamma <-|mmlb|manual", txt)
  cat(paste(trimws(txt[sort(unique(c(ix, ix+1, ix+2)))]), collapse=" | "), "\n")
}

# 2. empirical consequence: same seed, both arms, large N; compare the DEPENDENCE
#    (lag-1 within-cluster odds ratio and ICC on the latent scale proxy) and the
#    marginal prevalence. If the arms differed only in the tested beta, the
#    H0-vs-H1 difference in *within-cluster dependence* would be nil.
dep_stats <- function(d) {
  d <- d[order(d$id, d$time), ]
  s <- split(seq_len(nrow(d)), d$id)
  a <- b <- cc <- dd <- 0
  for (ix in s) {
    y <- d$binY[ix]; if (length(y) < 2) next
    for (j in 2:length(y)) {
      if (y[j-1]==1 && y[j]==1) a <- a+1
      if (y[j-1]==1 && y[j]==0) b <- b+1
      if (y[j-1]==0 && y[j]==1) cc <- cc+1
      if (y[j-1]==0 && y[j]==0) dd <- dd+1
    }
  }
  # cluster-level ICC of the binary outcome (ANOVA estimator)
  gm <- mean(d$binY)
  ni <- tapply(d$binY, d$id, length); mi <- tapply(d$binY, d$id, mean)
  msb <- sum(ni*(mi-gm)^2)/(length(ni)-1); msw <- sum((d$binY-mi[as.character(d$id)])^2)/(nrow(d)-length(ni))
  n0 <- (sum(ni)-sum(ni^2)/sum(ni))/(length(ni)-1)
  icc <- (msb-msw)/(msb+(n0-1)*msw)
  c(prev = gm, lag1_or = (a*dd)/max(b*cc,1), icc = icc,
    a=a,b=b,cc=cc,dd=dd, nobs = nrow(d))
}
set.seed(1)
res <- list(); k <- 0
for (m in 3:5) for (T in c(4, 10)) {
  for (lab in c("mmlb","manual")) {
    acc <- list()
    for (rep in 1:20) {
      set.seed(20260805L + 1000L*m + rep)   # production seeding scheme
      d <- make_all_data(600, T, rep, as.character(m), lab)
      acc[[rep]] <- d
    }
    D <- bind_rows(acc)
    k <- k+1; res[[k]] <- data.frame(model=m, ntime=T, arm=lab, t(dep_stats(D)))
  }
}
A <- bind_rows(res)
cat("\n--- GenBinaryY: marginal prevalence, lag-1 OR, ICC by arm (20 reps x N=600) ---\n")
print(as.data.frame(A %>% mutate(across(c(prev,lag1_or,icc), ~round(.x,4)))), row.names=FALSE)
cat("\n--- H1 minus H0 within the same (model, T) ---\n")
print(as.data.frame(A %>% select(model,ntime,arm,prev,lag1_or,icc) %>%
  tidyr::pivot_wider(names_from=arm, values_from=c(prev,lag1_or,icc)) %>%
  mutate(d_prev=round(prev_mmlb-prev_manual,4), d_lag1or=round(lag1_or_mmlb-lag1_or_manual,4),
         ratio_lag1or=round(lag1_or_mmlb/lag1_or_manual,4), d_icc=round(icc_mmlb-icc_manual,4)) %>%
  select(model,ntime,d_prev,d_lag1or,ratio_lag1or,d_icc)), row.names=FALSE)

cat("\n################ A2. the same check for GLMM-native and AR1 (D identical in both arms?)\n")
e_man <- new.env(); e_ar1 <- new.env()
sys.source(file.path(RD,"gen_glmm_data.R"), envir=e_man)
for (f in c("manual/generation_model3.R","manual/generation_model5.R")) sys.source(file.path(RD,f), envir=e_man)
cat("GLMM-native M3 D (both arms, hard-coded):\n"); print(eval(quote(D), envir=list2env(list(D=matrix(c(1,.01,.01,.001),2,byrow=TRUE)))))
cat("GLMM-native beta_c M3: mmlb", deparse(c(-2.32,0.116,0.232,0.232)), " manual", deparse(c(-2.32,0.116,0.232,0)), "\n")

cat("\n################ B/D. LRT vs Wald random-effects specification: how many cells?\n")
sp(file.path(SRC,"do_analysis.R"))
cat("do_analysis f_glmer (LRT):  binY ~ mean_rhs + (1|id)   [source line 68]\n")
for (m in 1:5) {
  a <- readLines(file.path(SRC, sprintf("analysis_model%d.R", m)))
  cat(sprintf("analysis_model%d glmer (Wald): %s\n", m, trimws(a[grep("binY ~", a)[length(grep("binY ~", a))]])))
}
NL <- c(10,15,20,seq(25,145,10)); TL <- c(2,3,4,6,8,10)
cat(sprintf("\ndesign points with a MISMATCHED LRT specification (M3-M5, all mechanisms): %d of %d (%.1f%%)\n",
  3*3*length(NL)*length(TL), 3*5*length(NL)*length(TL), 100*3/5))
cat(sprintf("(cell x hypothesis) affected: %d of %d\n", 2*3*3*length(NL)*length(TL), 2*3*5*length(NL)*length(TL)))

cat("\n################ C. GLMM-native: implied MARGINAL interaction by T, exact quadrature\n")
gh <- fastGHQuad::gaussHermiteData(80)
marg_p <- function(eta_c, V) sum(gh$w * plogis(eta_c + sqrt(2*V)*gh$x))/sqrt(pi)
implied <- function(bint, T, b0, bt, btr, V) {
  g <- expand.grid(time=0:(T-1), trt=0:1); X <- model.matrix(~time*trt, g)
  eta <- b0 + bt*g$time + btr*g$trt + bint*g$time*g$trt
  p <- vapply(eta, marg_p, numeric(1), V=V)
  beta <- rep(0, ncol(X))
  for (it in 1:300) { mu <- plogis(X%*%beta); W <- as.vector(mu*(1-mu))
    z <- X%*%beta + (p-mu)/pmax(W,1e-12)
    nb <- solve(crossprod(X, X*W), crossprod(X, W*z)); if (max(abs(nb-beta))<1e-13) {beta<-nb;break}; beta<-nb }
  as.vector(beta)[4]
}
g2 <- do.call(rbind, lapply(TL, function(T) data.frame(ntime=T,
  m_main = implied(0.232, T, -2.32, 0.116, 0.232, 1),
  m_m4   = implied(0.058, T, -2.32, 0.232, 0.232, 1))))
g2$dev_main <- g2$m_main - 0.20; g2$dev_m4 <- g2$m_m4 - 0.05
cat("random-intercept-only marginalisation (V = D[1,1] = 1):\n")
print(as.data.frame(g2 %>% mutate(across(everything(), ~round(.x,6)))), row.names=FALSE)
cat("\nNB: M3-M5 carry random SLOPES, so the observation-level latent variance is\n",
    "z'Dz, not D[1,1]. Implied marginal with the FULL D of M3 (var = 1 + 2*0.01*t + 0.001*t^2):\n")
implied_full <- function(bint, T, b0, bt, btr, Dm) {
  g <- expand.grid(time=0:(T-1), trt=0:1); X <- model.matrix(~time*trt, g)
  Z <- cbind(1, g$time); V <- rowSums((Z%*%Dm)*Z)
  eta <- b0 + bt*g$time + btr*g$trt + bint*g$time*g$trt
  p <- vapply(seq_along(eta), function(i) marg_p(eta[i], V[i]), numeric(1))
  beta <- rep(0, ncol(X))
  for (it in 1:300) { mu <- plogis(X%*%beta); W <- as.vector(mu*(1-mu))
    z <- X%*%beta + (p-mu)/pmax(W,1e-12)
    nb <- solve(crossprod(X, X*W), crossprod(X, W*z)); if (max(abs(nb-beta))<1e-13) {beta<-nb;break}; beta<-nb }
  as.vector(beta)[4]
}
Dm3 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
g3 <- do.call(rbind, lapply(TL, function(T) data.frame(ntime=T,
  m_M3 = implied_full(0.232, T, -2.32, 0.116, 0.232, Dm3))))
g3$dev <- g3$m_M3 - 0.20
print(as.data.frame(g3 %>% mutate(across(everything(), ~round(.x,6)))), row.names=FALSE)

cat("\n################ C2. AR1 conditional truth: SPEC used in truth_ar1_conditional.R vs the\n")
cat("random-effects structure the AR1 generator ACTUALLY used\n")
a <- readLines(file.path(RD,"ar1/generation_models_ar1.R"))
cat("generator, M3-M5 re.formula / mean.formula:\n")
print(grep("\\.ar1_generate\\(", a, value=TRUE))
b <- readLines(file.path(SRC,"postrun/truth_ar1_conditional.R"))
cat("\ntruth_ar1_conditional.R SPEC:\n"); print(grep("re = ~", b, value=TRUE))

sys.source(file.path(RD,"gen_ar1_data.R"), envir=globalenv())
sys.source(file.path(RD,"ar1/generation_models_ar1.R"), envir=globalenv())
mk_design <- function(N,T,ex) { set.seed(20260805L+1000L*3+1); .ar1_design(N,T,ex) }
cat("\nrecomputed conditional time:trt with the ACTUAL re.formula (~1+time) vs the\n",
    "SPEC used in F_truth_ar1.csv (~1+age), M3, N=75:\n")
out <- list(); k<-0
for (T in TL) {
  d <- mk_design(75,T,"age")
  D3 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
  bc_true <- calibrate_beta_c_ar1(~time*trt, ~1+time, .ar1_target(3,"mmlb"), D3, AR1_SIGMA2, d)
  bc_spec <- calibrate_beta_c_ar1(~time*trt, ~1+age,  .ar1_target(3,"mmlb"), D3, AR1_SIGMA2, d)
  k<-k+1; out[[k]] <- data.frame(ntime=T, cond_actual=as.numeric(bc_true)[4],
                                 cond_spec_in_Ftruth=as.numeric(bc_spec)[4], hardcoded=0.232)
}
O <- bind_rows(out); O$pct_err_hardcode <- 100*(O$hardcoded-O$cond_actual)/O$cond_actual
print(as.data.frame(O %>% mutate(across(everything(), ~round(.x,5)))), row.names=FALSE)
cat("\nsame for M1 (where the two SPECs agree, ~1):\n")
o1 <- do.call(rbind, lapply(TL, function(T) { set.seed(20260805L+1000L+1); d <- .ar1_design(75,T,"none")
  bc <- calibrate_beta_c_ar1(~time*trt, ~1, .ar1_target(1,"mmlb"), 1.0, AR1_SIGMA2, d)
  data.frame(ntime=T, cond=as.numeric(bc)[4]) }))
o1$pct_err <- 100*(0.232-o1$cond)/o1$cond
print(as.data.frame(o1 %>% mutate(across(everything(), ~round(.x,5)))), row.names=FALSE)
cat("\nAR1 H0 arm: is the calibrated beta_c identical to H1 except in the tested term?\n")
d <- mk_design(75,4,"age"); D3 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
b1 <- as.numeric(calibrate_beta_c_ar1(~time*trt, ~1+time, .ar1_target(3,"mmlb"), D3, AR1_SIGMA2, d))
b0 <- as.numeric(calibrate_beta_c_ar1(~time*trt, ~1+time, .ar1_target(3,"manual"), D3, AR1_SIGMA2, d))
print(rbind(H1=round(b1,6), H0=round(b0,6), diff=round(b1-b0,6)))

cat("\n################ SEEDS: is the seed a function of N, T, label, mechanism?\n")
cat("run_full.R line 118: set.seed(20260805L + 1000L * model + i)  -- i = replicate index\n")
cat("=> independent of N, ntime, label and dgp. Consequences, measured:\n")
sim_design <- function(N,T,m,lab,i) { set.seed(20260805L+1000L*m+i)
  d <- make_all_data(N,T,i,as.character(m),lab); d }
d10 <- sim_design(10,4,3,"mmlb",1); d15 <- sim_design(15,4,3,"mmlb",1)
cat(sprintf("N=10 vs N=15, m3 T4 rep1: trt of first 10 subjects identical? %s\n",
  identical(tapply(d10$trt,d10$id,`[`,1)[1:10], tapply(d15$trt,d15$id,`[`,1)[1:10])))
cat(sprintf("   cluster sizes of first 10 subjects identical? %s\n",
  identical(as.vector(table(d10$id))[1:10], as.vector(table(d15$id))[1:10])))
dH1 <- sim_design(45,4,3,"mmlb",7); dH0 <- sim_design(45,4,3,"manual",7)
cat(sprintf("H1 vs H0 same cell: design (id,time,trt) identical? %s ; outcome identical? %s\n",
  identical(dH1[,c("id","time","trt")], dH0[,c("id","time","trt")]), identical(dH1$binY, dH0$binY)))
cat(sprintf("seed range per model: model m uses [%d+1000m, %d+1000m+999]; overlap across models: %s\n",
  20260805, 20260805, "none (disjoint by construction, 1000 replicates per model)"))
cat("collision across (N,T,label,dgp) for the same (model, i): TOTAL -- same seed reused\n")
cat(sprintf("  number of cells sharing each seed: 16 N x 6 T x 2 labels x 3 dgp = %d\n", 16*6*2*3))

cat("\n################ TRUNCATION: 30%% of clusters lose follow-up\n")
tr <- do.call(rbind, lapply(TL, function(T) {
  ns <- c()
  for (i in 1:200) { set.seed(20260805L+1000L*3+i)
    N <- 500; nclust <- ifelse(runif(N) < 0.7, T, sample(seq_len(T-1), N, replace=TRUE)); ns <- c(ns, nclust) }
  data.frame(ntime=T, mean_obs=mean(ns), pct_full=100*mean(ns==T), pct_single=100*mean(ns==1),
             min_obs=min(ns)) }))
print(as.data.frame(tr %>% mutate(across(everything(), ~round(.x,3)))), row.names=FALSE)
cat("\nDONE_A\n")
