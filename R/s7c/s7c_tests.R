# =============================================================================
# s7c_tests.R -- S7c candidate tests for H0: time:trt = 0 (logit, marginal), with intercept, time, trt free.
# A flip       sign-flip of per-subject effective score contributions (Hemerik, Goeman, Finos 2020)
# B subj_t     two-stage subject-level analysis: residual time-trend r_i from the reduced model, Student t, N - 2 df
# C subj_perm  permutation of arm labels among r_i, Welch statistic (studentized permutation)
# D naive_perm permutation of the subject-level treatment, refit of the full model, |z| with CR0 sandwich (comparator)
# E naive_slope per-subject OLS slope of y on time, Student t between arms (comparator, probability scale)
# All resampling uses the RNG state set by the caller (set.seed(20261022 + i) after data generation).
# =============================================================================
B_FLIP <- 1999L; B_PERM <- 1999L; B_NAIVE <- 999L; ENUM_FLIP_MAXN <- 12L
COEF_MAX <- 10

s7c_prep <- function(d) {
  d <- as.data.frame(d); d <- d[order(d$id, d$time), ]
  id <- as.integer(factor(d$id))
  list(id = id, y = as.numeric(d$binY), time = as.numeric(d$time), trt = as.numeric(d$trt), N = max(id),
       trt_s = as.numeric(tapply(d$trt, id, `[`, 1)), n_i = as.integer(tabulate(id)))
}
# reduced (null) model: binY ~ time + trt, independence working correlation
null_fit <- function(P) {
  X0 <- cbind(1, P$time, P$trt)
  f <- tryCatch(suppressWarnings(glm.fit(X0, P$y, family = binomial(), control = glm.control(epsilon = 1e-12, maxit = 100))), error = function(e) NULL)
  if (is.null(f)) return(list(ok = FALSE, why = "null fit error"))
  b <- f$coefficients
  ok <- isTRUE(f$converged) && all(is.finite(b)) && max(abs(b)) < COEF_MAX
  list(ok = ok, why = if (ok) NA_character_ else if (!isTRUE(f$converged)) "null fit not converged" else "null fit coefficient non-finite or >= 10",
       mu = f$fitted.values, beta = b, X0 = X0, converged = isTRUE(f$converged))
}
# per-subject raw and effective score contributions for time:trt, and residual summaries r_i
score_parts <- function(P, F) {
  x <- P$time * P$trt; e <- P$y - F$mu; w <- F$mu * (1 - F$mu)
  X0 <- F$X0
  I00 <- crossprod(X0, X0 * w); Ix0 <- colSums(X0 * (w * x))
  u_raw <- rowsum(x * e, P$id, reorder = TRUE)[, 1]
  u_nu  <- rowsum(X0 * e, P$id, reorder = TRUE)
  proj  <- solve(I00, Ix0)
  u_eff <- as.numeric(u_raw - u_nu %*% proj)
  r <- rowsum(P$time * e, P$id, reorder = TRUE)[, 1]
  list(u_raw = as.numeric(u_raw), u_eff = u_eff, r = as.numeric(r))
}
all_signs <- function(N) {   # 2^N x N matrix of +-1, first row all +1
  m <- as.matrix(expand.grid(rep(list(c(1, -1)), N))); dimnames(m) <- NULL; m
}
test_flip <- function(u) {
  N <- length(u); To <- abs(sum(u))
  if (!(sum(u^2) > 0)) return(list(p = NA_real_, why = "zero score variance", n_ref = 0L))
  if (N <= ENUM_FLIP_MAXN) { G <- all_signs(N); Ts <- abs(as.numeric(G %*% u)); p <- mean(Ts >= To * (1 - 1e-10)); nref <- nrow(G) }
  else { G <- matrix(sample(c(-1, 1), B_FLIP * N, replace = TRUE), B_FLIP, N); Ts <- abs(as.numeric(G %*% u))
         p <- (1 + sum(Ts >= To * (1 - 1e-10))) / (B_FLIP + 1); nref <- B_FLIP + 1L }
  list(p = p, stat = sum(u), why = NA_character_, n_ref = as.integer(nref))
}
student_t <- function(r, g) {   # g in {0,1}
  r1 <- r[g == 1]; r0 <- r[g == 0]; n1 <- length(r1); n0 <- length(r0)
  if (n1 < 2 || n0 < 2) return(list(t = NA_real_, df = NA_real_, why = "fewer than 2 subjects in an arm"))
  sp2 <- ((n1 - 1) * var(r1) + (n0 - 1) * var(r0)) / (n1 + n0 - 2)
  if (!(sp2 > 0)) return(list(t = NA_real_, df = NA_real_, why = "zero pooled variance"))
  list(t = (mean(r1) - mean(r0)) / sqrt(sp2 * (1 / n1 + 1 / n0)), df = n1 + n0 - 2, why = NA_character_)
}
welch_mat <- function(r, G) {   # G: K x N 0/1 assignment matrix; Welch t for every row
  n1 <- rowSums(G); n0 <- ncol(G) - n1
  s1 <- as.numeric(G %*% r); q1 <- as.numeric(G %*% r^2); s0 <- sum(r) - s1; q0 <- sum(r^2) - q1
  m1 <- s1 / n1; m0 <- s0 / n0
  v1 <- (q1 - n1 * m1^2) / (n1 - 1); v0 <- (q0 - n0 * m0^2) / (n0 - 1)
  v1[v1 < 0] <- 0; v0[v0 < 0] <- 0
  se <- sqrt(v1 / n1 + v0 / n0); d <- m1 - m0
  out <- d / se; out[se == 0 & abs(d) < 1e-14] <- 0; out[se == 0 & abs(d) >= 1e-14] <- Inf * sign(d[se == 0 & abs(d) >= 1e-14])
  out
}
test_subj_perm <- function(r, g) {
  N <- length(r); n1 <- sum(g)
  if (n1 < 2 || N - n1 < 2) return(list(p = NA_real_, why = "fewer than 2 subjects in an arm", n_ref = 0L))
  if (!(var(r) > 0)) return(list(p = NA_real_, why = "zero variance of r", n_ref = 0L))
  to <- abs(welch_mat(r, matrix(g, 1)))
  if (choose(N, n1) <= B_PERM) {
    cb <- combn(N, n1); G <- matrix(0, ncol(cb), N); G[cbind(rep(seq_len(ncol(cb)), each = n1), as.vector(cb))] <- 1
    ts <- abs(welch_mat(r, G)); p <- mean(ts >= to * (1 - 1e-10)); nref <- nrow(G)
  } else {
    G <- t(vapply(seq_len(B_PERM), function(b) { z <- numeric(N); z[sample.int(N, n1)] <- 1; z }, numeric(N)))
    ts <- abs(welch_mat(r, G)); p <- (1 + sum(ts >= to * (1 - 1e-10))) / (B_PERM + 1); nref <- B_PERM + 1L
  }
  list(p = p, stat = welch_mat(r, matrix(g, 1)), why = NA_character_, n_ref = as.integer(nref))
}
# full model, independence working correlation, CR0 sandwich z for time:trt
z_full <- function(y, time, trt, id) {
  X <- cbind(1, time, trt, time * trt)
  f <- tryCatch(suppressWarnings(glm.fit(X, y, family = binomial())), error = function(e) NULL)
  if (is.null(f) || any(!is.finite(f$coefficients))) return(NA_real_)
  mu <- f$fitted.values; w <- mu * (1 - mu)
  Bi <- tryCatch(solve(crossprod(X, X * w)), error = function(e) NULL); if (is.null(Bi)) return(NA_real_)
  U <- rowsum(X * (y - mu), id); V <- Bi %*% crossprod(U) %*% Bi
  if (!(V[4, 4] > 0)) return(NA_real_)
  f$coefficients[4] / sqrt(V[4, 4])
}
test_naive_perm <- function(P) {
  zo <- z_full(P$y, P$time, P$trt, P$id)
  if (!is.finite(zo)) return(list(p = NA_real_, why = "observed z non-finite", n_ref = 0L, n_bad = NA_integer_))
  zs <- vapply(seq_len(B_NAIVE), function(b) { ts <- sample(P$trt_s); z_full(P$y, P$time, ts[P$id], P$id) }, 0)
  ok <- is.finite(zs)
  list(p = (1 + sum(abs(zs[ok]) >= abs(zo) * (1 - 1e-10))) / (sum(ok) + 1), stat = zo, why = NA_character_,
       n_ref = as.integer(sum(ok) + 1L), n_bad = as.integer(sum(!ok)))
}
test_naive_slope <- function(P) {
  sl <- vapply(seq_len(P$N), function(i) { k <- P$id == i; if (sum(k) < 2) return(NA_real_)
    tt <- P$time[k]; yy <- P$y[k]; sum((tt - mean(tt)) * (yy - mean(yy))) / sum((tt - mean(tt))^2) }, 0)
  has <- is.finite(sl); st <- student_t(sl[has], P$trt_s[has])
  p <- if (is.finite(st$t)) 2 * pt(-abs(st$t), st$df) else NA_real_
  list(p = p, stat = st$t, why = st$why, n_noslope = sum(!has),
       diff = if (is.finite(st$t)) mean(sl[has & P$trt_s == 1]) - mean(sl[has & P$trt_s == 0]) else NA_real_)
}
# all five tests on one dataset; returns a one-row data.frame per test (long format)
run_tests <- function(d, timing = FALSE) {
  P <- s7c_prep(d); F <- null_fit(P); g <- P$trt_s
  tm <- c(); tic <- function() proc.time()[["elapsed"]]
  row <- function(test, p, stat = NA_real_, why = NA_character_, n_ref = NA_integer_, extra = list())
    data.frame(test = test, p = p, stat = stat, why = why, n_ref = n_ref,
               n_noslope = if (is.null(extra$n_noslope)) NA_integer_ else extra$n_noslope,
               slope_diff = if (is.null(extra$diff)) NA_real_ else extra$diff,
               n_bad = if (is.null(extra$n_bad)) NA_integer_ else extra$n_bad, stringsAsFactors = FALSE)
  out <- list()
  t0 <- tic()
  if (F$ok) {
    S <- score_parts(P, F)
    a <- test_flip(S$u_eff); t1 <- tic()
    b <- student_t(S$r, g); pb <- if (is.finite(b$t)) 2 * pt(-abs(b$t), b$df) else NA_real_; t2 <- tic()
    cc <- test_subj_perm(S$r, g); t3 <- tic()
    out$A <- row("flip", a$p, if (is.null(a$stat)) NA_real_ else a$stat, a$why, a$n_ref)
    out$B <- row("subj_t", pb, b$t, b$why)
    out$C <- row("subj_perm", cc$p, if (is.null(cc$stat)) NA_real_ else cc$stat, cc$why, cc$n_ref)
    tm <- c(null_and_flip = t1 - t0, subj_t = t2 - t1, subj_perm = t3 - t2)
  } else {
    for (k in c("flip", "subj_t", "subj_perm")) out[[k]] <- row(k, NA_real_, why = F$why)
    t3 <- tic()
  }
  dd <- test_naive_perm(P); t4 <- tic()
  ee <- test_naive_slope(P); t5 <- tic()
  out$D <- row("naive_perm", dd$p, if (is.null(dd$stat)) NA_real_ else dd$stat, dd$why, dd$n_ref, list(n_bad = dd$n_bad))
  out$E <- row("naive_slope", ee$p, ee$stat, ee$why, extra = list(n_noslope = ee$n_noslope, diff = ee$diff))
  res <- do.call(rbind, out); rownames(res) <- NULL
  res$null_ok <- F$ok; res$N_subj <- P$N; res$n1 <- sum(g)
  if (timing) attr(res, "timing") <- c(tm, naive_perm = t4 - t3, naive_slope = t5 - t4)
  res
}
