# s7d_check_formulas.R -- independent hand-coded check of MD / KC on the geeglm fit (S7d preparation).
# KC per Kauermann and Carroll (2001): residuals premultiplied by the principal inverse square root of (I - H_ii),
# H_ii = D_i B^{-1} D_i' V_i^{-1} (non-symmetric). Principal root computed as Q diag(l^{-1/2}) Q^{-1} from eigen(I - H).
# geesmv's variant: same eigen() but Q' instead of Q^{-1} (valid only for symmetric matrices).
source(path.expand("~/MMM/R/s3/s3_common.R"))
.libPaths(c(.libPaths(), path.expand("~/R/s7d_probe_lib")))
suppressMessages({ library(clubSandwich); library(fewclust) })
env <- load_env("GenBinaryY", RUNDIR_V4)
hand <- function(fit, d) {
  X <- model.matrix(fit); mu <- fit$fitted.values; a <- mu * (1 - mu); y <- fit$y
  al <- unname(fit$geese$alpha); cs <- fit$corstr; ids <- split(seq_along(y), d$id)
  Rm <- function(n) switch(cs, independence = diag(n), exchangeable = { m <- matrix(al, n, n); diag(m) <- 1; m }, ar1 = al^abs(outer(1:n, 1:n, "-")))
  cl <- lapply(ids, function(ix) { n <- length(ix); Di <- X[ix, , drop = FALSE] * a[ix]; s <- sqrt(a[ix]); Vi <- (s %o% s) * Rm(n)
    list(D = Di, Vinv = solve(Vi), e = y[ix] - mu[ix]) })
  B <- Reduce(`+`, lapply(cl, function(z) t(z$D) %*% z$Vinv %*% z$D)); Bi <- solve(B)
  meat <- function(f) { M <- Reduce(`+`, lapply(cl, function(z) { u <- t(z$D) %*% z$Vinv %*% f(z); u %*% t(u) })); Bi %*% M %*% Bi }
  IH <- function(z) diag(nrow(z$D)) - z$D %*% Bi %*% t(z$D) %*% z$Vinv
  pr <- function(A) { E <- eigen(A); Re(E$vectors %*% diag(1 / sqrt(Re(E$values)), length(E$values)) %*% solve(E$vectors)) }
  gs <- function(A) { E <- eigen(A); d <- Re(E$values); d <- (d + abs(d)) / 2; d2 <- 1 / sqrt(d); d2[d == 0] <- 0; Re(E$vectors %*% diag(d2, length(d)) %*% t(E$vectors)) }
  k <- match("time:trt", colnames(X))
  c(md = sqrt(meat(function(z) solve(IH(z), z$e))[k, k]), kc_principal = sqrt(meat(function(z) pr(IH(z)) %*% z$e)[k, k]),
    kc_geesmv_formula = sqrt(meat(function(z) gs(IH(z)) %*% z$e)[k, k]))
}
out <- list()
for (cs in c("exchangeable", "ar1", "independence")) for (N in c(10L, 25L)) for (TT in c(3L, 10L)) for (i in 1:4) {
  set_rep_seed(1L, i); d <- as.data.frame(gen_data(env, "M1", 1L, N, TT, i, "manual")); d <- d[order(d$id, d$time), ]
  fit <- tryCatch(geepack::geeglm(binY ~ time * trt, id = id, data = d, family = binomial, corstr = cs, scale.fix = TRUE), error = function(e) NULL)
  if (is.null(fit) || max(abs(coef(fit))) > 10) next
  k <- match("time:trt", names(coef(fit)))
  h <- tryCatch(hand(fit, d), error = function(e) c(md = NA, kc_principal = NA, kc_geesmv_formula = NA))
  f <- sapply(c("md", "kc"), function(ty) tryCatch(sqrt(fewclust::gee_vcov(fit, ty)[k, k]), error = function(e) NA))
  cr <- sapply(c("CR2", "CR3"), function(ty) tryCatch(sqrt(as.matrix(clubSandwich::vcovCR(fit, cluster = d$id, type = ty))[k, k]), error = function(e) NA))
  sm <- if (cs == "ar1") NA else tryCatch(sqrt(geesmv::GEE.var.kc(binY ~ time * trt, id = "id", family = binomial, data = d, corstr = cs)$cov.beta[k]), error = function(e) NA)
  out[[length(out) + 1]] <- data.frame(cs, N, TT, i, t(h), fc_md = f[["md"]], fc_kc = f[["kc"]], CR2 = cr[["CR2"]], CR3 = cr[["CR3"]], geesmv_kc = sm)
}
R <- do.call(rbind, out); print(signif(R[, -1][, -(1:3)], 6)); print(R$cs)
dir.create(path.expand("~/MMM/results_v4/S7d/probe"), showWarnings = FALSE, recursive = TRUE)
write.csv(R, path.expand("~/MMM/results_v4/S7d/probe/check_formulas.csv"), row.names = FALSE)
