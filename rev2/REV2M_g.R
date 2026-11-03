# REV2M_g.R -- G: hand recomputation of MD / KC / FG from the estimating-equation
# building blocks on real fits, compared against geesmv on the SAME fit.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(geepack); library(geesmv); library(MMLB)})
SRC <- path.expand("~/MMM/R")
for (m in c(1,3)) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=globalenv())
sys.source(file.path(SRC,"generate_data.R"), envir=globalenv())
sys.source(file.path(SRC,"make_all_data.R"), envir=globalenv())

# ---- building blocks -------------------------------------------------------
blocks <- function(fit, d, formula) {
  X <- model.matrix(formula, d); y <- d$binY; idv <- d$id
  b <- coef(fit); p <- length(b)
  mu <- as.vector(plogis(X %*% b)); a <- mu*(1-mu)
  alpha <- unname(summary(fit)$corr[1,1])
  us <- unique(idv)
  I <- matrix(0,p,p); Mi <- vector("list", length(us)); Di <- Vi <- ei <- Mi
  for (k in seq_along(us)) {
    ix <- which(idv==us[k]); n_i <- length(ix)
    Xi <- X[ix,,drop=FALSE]; ai <- a[ix]
    Dk <- diag(ai, n_i) %*% Xi                     # d mu_i / d beta
    R  <- matrix(alpha, n_i, n_i); diag(R) <- 1
    A2 <- diag(sqrt(ai), n_i)
    Vk <- A2 %*% R %*% A2
    Vinv <- solve(Vk)
    Mk <- t(Dk) %*% Vinv %*% Dk
    I <- I + Mk
    Di[[k]] <- Dk; Vi[[k]] <- Vinv; ei[[k]] <- y[ix]-mu[ix]; Mi[[k]] <- Mk
  }
  Iinv <- solve(I)
  list(Iinv=Iinv, Di=Di, Vi=Vi, ei=ei, Mi=Mi, p=p, us=us, X=X, alpha=alpha)
}
msqrt_inv <- function(M) {                 # symmetric inverse square root
  e <- eigen((M+t(M))/2, symmetric=TRUE)
  e$vectors %*% diag(1/sqrt(pmax(e$values,1e-12))) %*% t(e$vectors)
}
sandwiches <- function(B) {
  p <- B$p; Iinv <- B$Iinv
  mLZ <- mMD <- mKC <- mFG <- matrix(0,p,p)
  for (k in seq_along(B$Di)) {
    Dk <- B$Di[[k]]; Vinv <- B$Vi[[k]]; ek <- B$ei[[k]]
    n_i <- length(ek)
    g <- t(Dk) %*% Vinv %*% ek
    mLZ <- mLZ + g %*% t(g)
    H <- Dk %*% Iinv %*% t(Dk) %*% Vinv            # observation-level leverage
    Im <- diag(n_i) - H
    # Mancl-DeRouen: (I - H)^{-1} e
    eMD <- tryCatch(solve(Im, ek), error=function(e) rep(NA,n_i))
    gMD <- t(Dk) %*% Vinv %*% eMD; mMD <- mMD + gMD %*% t(gMD)
    # Kauermann-Carroll: (I - H)^{-1/2} e  (symmetric root of a NON-symmetric Im)
    eKC <- tryCatch(as.vector(msqrt_inv(Im) %*% ek), error=function(e) rep(NA,n_i))
    gKC <- t(Dk) %*% Vinv %*% eKC; mKC <- mKC + gKC %*% t(gKC)
    # Fay-Graubard: diagonal adjustment in PARAMETER space, bounded by r
    r <- 0.75
    dd <- diag(Iinv %*% B$Mi[[k]])
    Hfg <- diag(1/sqrt(pmax(1 - pmin(r, dd), 1e-12)), p)
    gFG <- Hfg %*% g; mFG <- mFG + gFG %*% t(gFG)
  }
  list(LZ = Iinv %*% mLZ %*% Iinv, MD = Iinv %*% mMD %*% Iinv,
       KC = Iinv %*% mKC %*% Iinv, FG = Iinv %*% mFG %*% Iinv)
}

one <- function(N, T, model=1, rep=1) {
  set.seed(20260805L + 1000L*model + rep)
  d <- make_all_data(N, T, rep, as.character(model), "mmlb")
  d <- as.data.frame(d); d <- d[order(d$id),]
  f <- binY ~ time*trt
  fit <- tryCatch(geeglm(f, id=factor(id), data=d, family=binomial,
                         corstr="exchangeable", scale.fix=TRUE), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  k <- match("time:trt", names(coef(fit)))
  B <- blocks(fit, d, f); S <- sandwiches(B)
  pk <- lapply(c(md="md", kc="kc", fg="fg"), function(nm)
    tryCatch(do.call(getExportedValue("geesmv", paste0("GEE.var.",nm)),
             list(f, id="id", family=binomial, data=d, corstr="exchangeable")),
             error=function(e) NULL))
  se_pkg <- sapply(pk, function(v) if (is.null(v)) NA else sqrt(v$cov.beta[k]))
  se_geepack <- summary(fit)$coefficients[k,2]
  data.frame(N=N, ntime=T, model=model, rep=rep, n_obs=nrow(d), n_events=sum(d$binY),
             alpha=B$alpha,
             se_lz_hand=sqrt(S$LZ[k,k]), se_lz_geepack=se_geepack,
             se_md_hand=sqrt(S$MD[k,k]), se_md_pkg=se_pkg["md"],
             se_kc_hand=sqrt(S$KC[k,k]), se_kc_pkg=se_pkg["kc"],
             se_fg_hand=sqrt(S$FG[k,k]), se_fg_pkg=se_pkg["fg"])
}

cat("===== single fit, N=15 T=4, GenBinaryY M1, replicate 1 =====\n")
r0 <- one(15,4,1,1); print(as.data.frame(t(round(unlist(r0[,-(1:4)]),6))))
cat(sprintf("\nhand/package ratio  MD %.6f  KC %.6f  FG %.6f\n",
  r0$se_md_hand/r0$se_md_pkg, r0$se_kc_hand/r0$se_kc_pkg, r0$se_fg_hand/r0$se_fg_pkg))
cat(sprintf("hand/uncorrected    MD %.6f  KC %.6f  FG %.6f  (LZ hand/geepack %.6f)\n",
  r0$se_md_hand/r0$se_lz_hand, r0$se_kc_hand/r0$se_lz_hand, r0$se_fg_hand/r0$se_lz_hand,
  r0$se_lz_hand/r0$se_lz_geepack))

cat("\n===== grid of cluster counts, T=4, M1, 3 replicates each =====\n")
G <- expand.grid(N=c(10,15,20,25,45,65,85,105,145), rep=1:3)
res <- bind_rows(lapply(seq_len(nrow(G)), function(i) tryCatch(one(G$N[i],4,1,G$rep[i]), error=function(e) NULL)))
res <- res %>% mutate(
  r_md_pkg = se_md_hand/se_md_pkg, r_kc_pkg = se_kc_hand/se_kc_pkg, r_fg_pkg = se_fg_hand/se_fg_pkg,
  r_md_unc = se_md_hand/se_lz_hand, r_kc_unc = se_kc_hand/se_lz_hand, r_fg_unc = se_fg_hand/se_lz_hand,
  r_kcpkg_unc = se_kc_pkg/se_lz_geepack, r_mdpkg_unc = se_md_pkg/se_lz_geepack,
  r_fgpkg_unc = se_fg_pkg/se_lz_geepack)
print(as.data.frame(res %>% select(N,rep,n_events,r_md_pkg,r_kc_pkg,r_fg_pkg) %>%
  mutate(across(where(is.numeric), ~round(.x,5)))), row.names=FALSE)
cat("\n--- max |hand/package - 1| over the grid ---\n")
cat(sprintf("MD %.5f | KC %.5f | FG %.5f\n",
  max(abs(res$r_md_pkg-1),na.rm=TRUE), max(abs(res$r_kc_pkg-1),na.rm=TRUE),
  max(abs(res$r_fg_pkg-1),na.rm=TRUE)))
cat("\n--- ratio to the UNCORRECTED sandwich (hand), by N ---\n")
print(as.data.frame(res %>% group_by(N) %>% summarise(
  md=round(mean(r_md_unc),4), kc=round(mean(r_kc_unc),4), fg=round(mean(r_fg_unc),4), .groups="drop")),
  row.names=FALSE)
cat("\n--- ratio to the UNCORRECTED sandwich (PACKAGE), by N ---\n")
print(as.data.frame(res %>% group_by(N) %>% summarise(
  md=round(mean(r_mdpkg_unc),4), kc=round(mean(r_kcpkg_unc),4), fg=round(mean(r_fgpkg_unc),4),
  .groups="drop")), row.names=FALSE)
cat(sprintf("\nFG hand/uncorrected range: %.4f - %.4f | KC hand: %.4f - %.4f | KC pkg: %.4f - %.4f\n",
  min(res$r_fg_unc,na.rm=TRUE), max(res$r_fg_unc,na.rm=TRUE),
  min(res$r_kc_unc,na.rm=TRUE), max(res$r_kc_unc,na.rm=TRUE),
  min(res$r_kcpkg_unc,na.rm=TRUE), max(res$r_kcpkg_unc,na.rm=TRUE)))
cat(sprintf("max package KC/uncorrected ratio (looking for the 10^21 claim): %.6g\n",
  max(res$r_kcpkg_unc, na.rm=TRUE)))
write.csv(res, path.expand("~/MMM/rev2/REV2M_hand_corrections.csv"), row.names=FALSE)
cat("\nDONE_G\n")
